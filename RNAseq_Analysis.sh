#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --ntasks=4
#SBATCH --mem=16g
#SBATCH --tmp=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=howex118@umn.edu

helpFunction()
{
   echo ""
   echo "Usage: $0 -f parameterFolder -i parameterInput -d parameterDNA -a parameterAnno"
   echo -e "\t-f Folder where RNA-sequencing analysis will take place"
   echo -e "\t-i Input directory containing paired-end zipped fastq files to be analyzed"
   echo -e "\t-d DNA reference file directory path"
   echo -e "\t-a Annotation file directory path"
   exit 1 # Exit script after printing help
}

while getopts "f:i:d:a:" opt
do
   case "$opt" in
      f ) parameterFolder="$OPTARG" ;;
      i ) parameterInput="$OPTARG" ;;
      d ) parameterDNA="$OPTARG" ;;
      a ) parameterAnno="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done
echo "$parameterFolder"
echo "$parameterInput"
echo "$parameterDNA"
echo "$parameterAnno"
# Print helpFunction in case parameters are empty
if [ -z "$parameterFolder" ] || [ -z "$parameterInput" ] || [ -z "$parameterDNA" ] || [ -z "$parameterAnno" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$parameterFolder"
echo "$parameterInput"
echo "$parameterDNA"
echo "$parameterAnno"

module load cutadapt/2.4
module load star/2.7.1a
module load fastqc/0.11.7

echo "Step 1: Transferring RNAseq files to analysis directory and running FastQC on them"

cp -t $parameterFolder $parameterInput/*.gz

cd $parameterFolder

mkdir trimmed_reads
mkdir Index
mkdir STAR
mkdir Genome_indices
mkdir Fastqc

echo "Step 2: Create sample list"
ls *_R1_001.fastq.gz | cut -f1 -d "." > samples_names_RNAseqR1.txt
ls *_R2_001.fastq.gz | cut -f1 -d "." > samples_names_RNAseqR2.txt


echo "Step 3: Running Fastqc"
paste samples_names_RNAseqR1.txt samples_names_RNAseqR2.txt | while read sampleR1 sampleR2; do

    fastqc -t 4 ${sampleR1}.fastq.gz -o $parameterFolder/Fastqc
    fastqc -t 4 ${sampleR2}.fastq.gz -o $parameterFolder/Fastqc

done


echo "Step 4: Trimming illumina adapters from files"
paste samples_names_RNAseqR1.txt samples_names_RNAseqR2.txt | while read sampleR1 sampleR2; do

    echo "On sample: $sampleR1 and $sampleR2 "

    cutadapt -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -m 30 \
    -u 1 \
    --quality-base=33 \
    --cores=8 \
    -o $parameterFolder/trimmed_reads/${sampleR1}_trimmed.fastq.gz -p $parameterFolder/trimmed_reads/${sampleR2}_trimmed.fastq.gz \
    $parameterFolder/${sampleR1}.fastq.gz $parameterFolder/${sampleR2}.fastq.gz \
    >> cutadapt_primer_trimming_stats.txt 2>&1

done

echo "Step 5: Align reads to Mycobacterium genome using STAR"
echo "Creating genome index"
STAR --runMode genomeGenerate \
     --genomeSAindexNbases 8 \
     --genomeDir $parameterFolder/Genome_indices \
     --runThreadN 8 \
     --genomeFastaFiles $parameterDNA \

echo "Aligning reads for each sample"

paste samples_names_RNAseqR1.txt samples_names_RNAseqR2.txt | while read sampleR1 sampleR2; do

    echo "On sample: $sampleR1 and $sampleR2 "

    STAR --genomeDir $parameterFolder/Genome_indices \
         --runThreadN 8 \
         --readFilesIn $parameterFolder/trimmed_reads/${sampleR1}_trimmed.fastq.gz $parameterFolder/trimmed_reads/${sampleR2}_trimmed.fastq.gz \
         --readFilesCommand zcat \
         --alignIntronMax 1 \
         --limitBAMsortRAM 1172893133 \
         --outFileNamePrefix $parameterFolder/STAR/${sampleR1}_trimmed_STAR_Aligned \
         --outSAMtype BAM SortedByCoordinate

done

echo "Step 6: Counting reads with featureCounts"
echo "Begin counting"
/home/baughna/howex118/subread-2.0.3-source/bin/featureCounts -p --countReadPairs -s 2 -t gene -T 4 \
-g locus_tag \
--extraAttributes gene \
-a $parameterAnno \
-o RNAseq.featureCounts.txt \
$parameterFolder/STAR/*.bam
