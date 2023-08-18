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


module load cutadapt/2.4
module load star/2.7.1a
module load fastqc/0.11.7

echo "Step 1: Transferring RNAseq files to analysis directory and running FastQC on them"

cd $parameterFolder

mkdir Input_reads
mkdir trimmed_reads
mkdir STAR
mkdir Genome_indices
mkdir Fastqc


FastaRef=($parameterDNA/*.fasta)
AnnoRef=($parameterDNA/*.gtf)


echo "$FastaRef"
echo "$AnnoRef"


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

echo "Step 4: Trimming t-overhangs"
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

# Old BWA code, should work if needed to test something.

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