#!/bin/sh
#$ -N rnaseq
#$ -cwd
#$ -j y
#$ -o UGE-output/
#$ -m be
#$ -M username@nih.gov
#$ -pe threaded 6


#cd /hpcdata/{group}/{workingdir} 

# Copy input files
cp -r /hpcdata/scratch/rnaseq_lesson1 .
cd rnaseq_lesson1


# Which modules to use?  Load and unload as needed
# module load STAR/2.7.3a-goolf-1.7.20 
# module load samtools
# module load fastqc/0.11.8-Java-1.8.0_45
# module load multiqc/1.7-Python-3.5.5
# module load cutadapt/2.3-foss-2016b-Python-3.7.3
# module load picard/2.17.6
# module load HTSeq/0.9.1-goolf-1.7.20-Python-2.7.9

# Quality check

module load fastqc/0.11.8-Java-1.8.0_45
module load multiqc/1.7-Python-3.5.5

cd raw_data
fastqc *fastq
multiqc .

module unload fastqc/0.11.8-Java-1.8.0_45
module unload multiqc/1.7-Python-3.5.5

# Exercise 1: Trim first 10 bases, then rerun fastqc and examine output
module load cutadapt/2.3-foss-2016b-Python-3.7.3

mkdir trimmedreads

file="inputIDs.txt"
while read line
    do
    cat $line | cutadapt -u 10 -U 10 --minimum-length=25 -o trimmedreads/${line}.read1.trimmed.fastq -p trimmedreads/${line}.read2.trimmed.fastq ${line}_Build37-ErccTranscripts-chr22.read1.fastq ${line}_Build37-ErccTranscripts-chr22.read2.fastq
    echo "$line"
done <"$file"

module unload cutadapt/2.3-foss-2016b-Python-3.7.3
module load fastqc/0.11.8-Java-1.8.0_45
module load multiqc/1.7-Python-3.5.5

cd trimmedreads
fastqc *trimmed.fastq
multiqc -f .

module unload fastqc/0.11.8-Java-1.8.0_45
module unload multiqc/1.7-Python-3.5.5


# Create index and align
module load STAR/2.7.3a-goolf-1.7.20 

cd ../../ # you should now be in the directory `rnaseq_lesson1`
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir reference_data/chr22index \
--genomeFastaFiles reference_data/chr22.fa \
--sjdbGTFfile reference_data/genes.gtf \
--sjdbOverhang 89 \
--genomeSAindexNbases 11


# Exercise 2: Align reads
# do HBR_Rep1
# STAR --genomeDir ../reference_data/chr22index \
# --runThreadN 6 \
# --readFilesCommand trimmedreads/HBR_Rep1_ERCC-Mix2.read1.trimmed.fastq HBR_Rep1_ERCC-Mix2.read2.trimmed.fastq \
# --outFileNamePrefix results/HBR_Rep1_test \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes Standard 

# # or UHR_Rep3
# STAR --genomeDir ../reference_data/chr22index \
# --runThreadN 6 \
# --readFilesCommand trimmedreads/UHR_Rep3_ERCC-Mix1.read1.trimmed.fastq trimmedreads/UHR_Rep3_ERCC-Mix1.read2.trimmed.fastq \
# --outFileNamePrefix results/UHR_Rep3_test \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes Standard

# or Use a Loop to align all
cd raw_data
file="inputIDs.txt"
while read line
    do
    STAR \
    --outSAMtype BAM SortedByCoordinate \
    --genomeDir ../reference_data/chr22index \
    --readFilesIn trimmedreads/${line}.read1.trimmed.fastq trimmedreads/${line}.read2.trimmed.fastq \
    --runThreadN 6 \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFileNamePrefix ../results/${line}_
done <"$file"

module unload STAR/2.7.3a-goolf-1.7.20 

# rename files
cd ../results
cp HBR_Rep1*bam HBR_Rep1.bam
cp HBR_Rep2*bam HBR_Rep2.bam
cp HBR_Rep3*bam HBR_Rep3.bam
cp UHR_Rep1*bam UHR_Rep1.bam
cp UHR_Rep2*bam UHR_Rep2.bam
cp UHR_Rep3*bam UHR_Rep3.bam


#**Exercise 3: Prepare files to view in IGV**

# If you have not yet aligned all 6 samples, copy the directory rnaseq_bams_lesson1 to the working directory
#cp -r /hpcdata/scratch/rnaseq_bams_lesson1 .

module load picard/2.17.6

#cd rnaseq_bams_lesson1
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar MergeSamFiles OUTPUT=UHR.bam INPUT=UHR_Rep1.bam INPUT=UHR_Rep2.bam INPUT=UHR_Rep3.bam
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar MergeSamFiles OUTPUT=HBR.bam INPUT=HBR_Rep1.bam INPUT=HBR_Rep2.bam INPUT=HBR_Rep3.bam

module unload picard/2.17.6

# index the bam files
module load samtools

samtools index UHR.bam
samtools index HBR.bam

module unload samtools

#**Exercise 4: Create counts table and explore which genes have the most counts - preparation for differential abundance

module load HTSeq/0.9.1-goolf-1.7.20-Python-2.7.9


htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id UHR_Rep1.bam ../reference_data/genes.gtf > UHR_Rep1_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id UHR_Rep2.bam ../reference_data/genes.gtf > UHR_Rep2_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id UHR_Rep3.bam ../reference_data/genes.gtf > UHR_Rep3_gene.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id HBR_Rep1.bam ../reference_data/genes.gtf > HBR_Rep1_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id HBR_Rep2.bam ../reference_data/genes.gtf > HBR_Rep2_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id HBR_Rep3.bam ../reference_data/genes.gtf > HBR_Rep3_gene.tsv

join UHR_Rep1_gene.tsv UHR_Rep2_gene.tsv | join - UHR_Rep3_gene.tsv | join - HBR_Rep1_gene.tsv | join - HBR_Rep2_gene.tsv | join - HBR_Rep3_gene.tsv > gene_read_counts_table_all.tsv
echo "GeneID UHR_Rep1 UHR_Rep2 UHR_Rep3 HBR_Rep1 HBR_Rep2 HBR_Rep3" > header.txt
cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > gene_read_counts_table_all_final.tsv
rm -f gene_read_counts_table_all.tsv header.txt
head gene_read_counts_table_all_final.tsv

