#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --ntasks=40
#SBATCH --mem=32g
#SBATCH --tmp=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=howex118@umn.edu

helpFunction() {
   echo ""
   echo "Usage: $0 -f parameterFolder -i parameterInput -d parameterDNA -a parameterAnno"
   echo -e "\t-f Folder where RNA-sequencing analysis will take place"
   echo -e "\t-i Input directory containing paired-end zipped fastq files to be analyzed"
   echo -e "\t-d DNA reference files directory path"
   exit 1 # Exit script after printing help
}

while getopts "f:i:d:" ARGS; do
   case "$ARGS" in
   f) parameterFolder=$OPTARG ;;
   i) parameterInput=$OPTARG ;;
   d) parameterDNA=$OPTARG ;;
   ?) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterFolder" ] || [ -z "$parameterInput" ] || [ -z "$parameterDNA" ]; then
   echo "Some or all of the parameters are empty"
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$parameterFolder"
echo "$parameterInput"
echo "$parameterDNA"

module use ~/modulefiles.local
module load cutadapt/2.4
module load star/2.7.1a
module load fastqc/0.11.7
module load python3
# module load ribodetector
module load samtools
module load RSeQC
module load multiqc
module load R

echo "Step 1: Transferring RNAseq files to analysis directory and running FastQC on them"

cd $parameterFolder

mkdir Input_reads
mkdir trimmed_reads
mkdir STAR
mkdir Genome_indices
mkdir Fastqc
mkdir ribodepleted

FastaRef=($parameterDNA/*.fasta)
AnnoRef=($parameterDNA/*.gtf)
BedRef=($parameterDNA/*.bed)

echo "$FastaRef"
echo "$AnnoRef"
echo "$BedRef"

cp -t ${parameterFolder}/Input_reads ${parameterInput}/*.gz

(
   cd Input_reads

   echo "Step 2: Create sample list"
   ls *_R1_001.fastq.gz | cut -f1 -d "." >$parameterFolder/samples_names_RNAseqR1.txt
   ls *_R2_001.fastq.gz | cut -f1 -d "." >$parameterFolder/samples_names_RNAseqR2.txt
)
echo "Step 3: Running Fastqc"
paste samples_names_RNAseqR1.txt samples_names_RNAseqR2.txt | while read sampleR1 sampleR2; do

   fastqc -t 4 Input_reads/${sampleR1}.fastq.gz -o $parameterFolder/Fastqc
   fastqc -t 4 Input_reads/${sampleR2}.fastq.gz -o $parameterFolder/Fastqc

done

echo "Step 4: Trimming t-overhangs and removing rRNA from samples"
paste samples_names_RNAseqR1.txt samples_names_RNAseqR2.txt | while read sampleR1 sampleR2; do

   echo "Trimming onn sample: $sampleR1 and $sampleR2 "

   cutadapt -m 30 \
      -u 1 \
      -U 1 \
      --quality-base=33 \
      --cores=40 \
      -o $parameterFolder/trimmed_reads/${sampleR1}_trimmed.fastq.gz -p $parameterFolder/trimmed_reads/${sampleR2}_trimmed.fastq.gz \
      $parameterFolder/Input_reads/${sampleR1}.fastq.gz $parameterFolder/Input_reads/${sampleR2}.fastq.gz \
      >>cutadapt_primer_trimming_stats.txt 2>&1

done

# paste samples_names_RNAseqR1.txt samples_names_RNAseqR2.txt | while read sampleR1 sampleR2; do

#    echo "Removing rRNA in sample: $sampleR1 and $sampleR2 "

#    ribodetector_cpu -t 40 \
#       -l 50 \
#       -i $parameterFolder/trimmed_reads/${sampleR1}_trimmed.fastq.gz $parameterFolder/trimmed_reads/${sampleR2}_trimmed.fastq.gz \
#       -e norrna \
#       -r ribodepleted/${sampleR1}_rRNAonly.fastq ribodepleted/${sampleR2}_rRNAonly.fastq \
#       --chunk_size 256 \
#       -o ribodepleted/${sampleR1}_norRNA.fastq ribodepleted/${sampleR2}_norRNA.fastq
# done

echo "Step 5: Align reads to Mycobacterium genome using STAR"
echo "Creating genome index"
STAR --runMode genomeGenerate \
   --genomeSAindexNbases 8 \
   --genomeDir $parameterFolder/Genome_indices \
   --runThreadN 8 \
   --genomeFastaFiles $FastaRef

echo "Aligning reads for each sample"

paste samples_names_RNAseqR1.txt samples_names_RNAseqR2.txt | while read sampleR1 sampleR2; do

   echo "On sample: $sampleR1 and $sampleR2 "

   STAR --genomeDir $parameterFolder/Genome_indices \
      --runThreadN 8 \
      --readFilesIn $parameterFolder/trimmed_reads/${sampleR1}_trimmed.fastq.gz $parameterFolder/trimmed_reads/${sampleR2}_trimmed.fastq.gz \
      --alignIntronMax 1 \
      --limitBAMsortRAM 1172893133 \
      --outFileNamePrefix $parameterFolder/STAR/${sampleR1}_STAR \
      --outSAMtype BAM SortedByCoordinate

done

# echo "Step 5: Align reads to Mycobacterium genome using STAR"
# echo "Creating genome index"

# bwa index ${parameterFolder}/Genome_indices
# paste samples_names_RNAseqR1.txt samples_names_RNAseqR2.txt | while read sampleR1 sampleR2; do

#    echo "On sample: $sampleR1 and $sampleR2 "

#    bwa mem -t 40 $parameterFolder/Genome_indices/$parameterDNA $parameterFolder/ribodepeleted/${sampleR1}_norRNA.fastq $parameterFolder/ribodepeleted/${sampleR2}_norRNA.fastq | samtools sort -@40 -o $parameterFolder/STAR/${sampleR1}.bam -

# done

echo "Step 6: Counting reads with featureCounts"
echo "Begin counting"
/home/baughna/howex118/subread-2.0.3-source/bin/featureCounts -p --countReadPairs -s 2 -t gene -T 4 \
   -g locus_tag \
   --extraAttributes gene \
   -a $AnnoRef \
   -o RNAseq.featureCounts.txt \
   $parameterFolder/STAR/*.bam

echo "Step 7: QC Metrics"
echo "Create .bai indexes from .bam files for downstream analysis tools"

paste samples_names_RNAseqR1.txt | while read sampleR1; do

   echo "Indexing sample: $sampleR1.bam"
   samtools index $parameterFolder/STAR/${sampleR1}_STARAligned.sortedByCoord.out.bam

done

geneBody_coverage.py -l 500 -r $BedRef -i $parameterFolder/STAR/ -o genebodycoverageData

paste samples_names_RNAseqR1.txt | while read sampleR1; do

   echo "Completing functions on sample: $sampleR1"

   inner_distance.py -i $parameterFolder/STAR/${sampleR1}_STARAligned.sortedByCoord.out.bam -o ${sampleR1}_innerdistance -r $BedRef
   read_distribution.py -i $parameterFolder/STAR/${sampleR1}_STARAligned.sortedByCoord.out.bam -r $BedRef
   read_duplication.py -i $parameterFolder/STAR/${sampleR1}_STARAligned.sortedByCoord.out.bam -o ${sampleR1}_read_duplication

done

echo "Step 8: Run MultiQC"
multiqc $parameterFolder/ \
   $parameterFolder/STAR \
   # $parameterFolder/ribodepleted \
   $parameterFolder/Fastqc
