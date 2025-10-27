#!/bin/bash

# A script to perform the preprocessing of HUVEC ZEB1 KD RNA-seq data.
# This is largely based on our workflow publication: https://pubmed.ncbi.nlm.nih.gov/35099752/

# Create a main working directory and set the path to a variable
mkdir /~/cnv_paper/data/bulkrnaseq
main_path=/~/cnv_paper/data/bulkrnaseq

# Create directories to store the raw data as well as the outputs from processing.
# raw_data should be where the fastq files are stored.

mkdir $main_path/raw_data
mkdir $main_path/outputs
mkdir $main_path/outputs/fastqc
mkdir $main_path/outputs/cutadapt
mkdir $main_path/outputs/star
mkdir $main_path/outputs/featurecounts

# Specify variables for each directory for ease of use in commands.

data_dir=$main_path/raw_data
outputs_dir=$main_path/outputs
fastqc_dir=$main_path/outputs/fastqc
cutadapt_dir=$main_path/outputs/cutadapt
star_dir=$main_path/outputs/star
counts_dir=$main_path/outputs/featurecounts

# We assume that a genome has already been downloaded, indexed and saved in ~/star_index/human/starindex107HS/.
# If not, see our publication https://pubmed.ncbi.nlm.nih.gov/35099752/

# List of samples to condition below.
# The data is paired end, so there is a 1 and 2 for each sample.

# | Sample ID | Condition |
# | --------- | --------- |
# | NS2_1 | HUVEC_control\_1
# | NS2_2 | HUVEC_control\_1
# | NS3_1 | HUVEC_control\_2
# | NS3_2 | HUVEC_control\_2
# | NS8_1 | HUVEC_control\_3
# | NS8_2 | HUVEC_control\_3
# | ZEB1_4\_1 | HUVEC_ZEB1KD\_1
# | ZEB1_4\_2 | HUVEC_ZEB1KD\_1
# | ZEB1_5\_1 | HUVEC_ZEB1KD\_2
# | ZEB1_5\_2 | HUVEC_ZEB1KD\_2
# | ZEB1_6\_1 | HUVEC_ZEB1KD\_3
# | ZEB1_6\_2 | HUVEC_ZEB1KD\_3

samples="NS2_1 NS2_2 NS3_1 NS3_2 NS8_1 NS8_2 ZEB1_4_1 ZEB1_4_2 ZEB1_5_1 ZEB1_5_2 ZEB1_6_1 ZEB1_6_2"

# Note - we ran this analysis on unzipped fastq files.
# This analysis can be done on .fastq.gz files, but be aware of any necessary arguments.

# Make sure all the .fastq files are transferred into the raw_data directory.

cd $data_dir

# Run fastqc on all the files

fastqc -t 14 *.fastq -o $fastqc_dir

# Check the fastq files are properly formatted.with cutadapt

for i in *.fastq
do
    cutadapt -o /dev/null $i -j 14
done

# Complete a further check to make sure the FASTQ files are correctly paired-end, if the data is paired-end. Ignore if single-end.

for i in $samples
do
    cutadapt -o /dev/null -p /dev/null "$i"_1.fastq "$i"_2.fastq -j 14
done

# Run cutadapt. The parameters and adapters can be changed to suit based on fastqc results.
# Note this command is specifically for paired-end data. For single end see the cutadapt documentation.
# -j is the threads. Adjust accordingly.

for i in $samples
do
    cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -q 20 --minimum-length 1 -o $cutadapt_dir/"$i"_1_trimmed.fastq -p $cutadapt_dir/"$i"_2_trimmed.fastq "$i"_1.fastq "$i"_2.fastq -j 14
done

# Align with STAR to human genome.
# In our case, our indexed genome is in ~/star_index/human/starindex107HS/
# Adjust this path accordingly.
# --readFilesCommand zcat may be needed for .fastq.gz files - check this.

cd $cutadapt_dir

for i in $samples
do
    STAR --genomeDir ~/star_index/human/starindex107HS/ --readFilesCommand zcat --readFilesIn $cutadapt_dir/"$i"_1_trimmed.fq.gz $cutadapt_dir/"$i"_2_trimmed.fq.gz --runThreadN 14 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $star_dir/$i > $i_align.log
done

# At this point Qualimap can be run to check the strandedness of the data and analyse the alignment.
# This needs to be known for quantification, if it is not already known.
# See our publication for more information: https://pubmed.ncbi.nlm.nih.gov/35099752/

# Quantifying the aligned reads.

cd $star_dir

featureCounts -T 12 -p -s 0 -t exon -g gene_id -a ~/star_index/human/starindex107HS/Homo_sapiens.GRCh38.107.chr.gtf -o $counts_dir/HUVEC_ZEB1_KD.txt *.bam

# The output file from featurecounts can be loaded into DESeq2 in R for downstream analyses.
