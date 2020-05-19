# run this in an interactive session on Locus

qrsh
cp -r /hpcdata/scratch/rnaseq_bams_lesson1 .
cd rnaseq_bams_lesson1

# ls  # you should see:
# HBR_Rep1.bam  HBR_Rep2.bam  HBR_Rep3.bam  UHR_Rep1.bam  UHR_Rep2.bam  UHR_Rep3.bam

# Run Samtools flagstat
module load samtools

samtools flagstat HBR_Rep1.bam > HBR_Rep1_flagstat.txt
samtools flagstat HBR_Rep2.bam > HBR_Rep2_flagstat.txt
samtools flagstat HBR_Rep3.bam > HBR_Rep3_flagstat.txt

samtools flagstat UHR_Rep1.bam > UHR_Rep1_flagstat.txt
samtools flagstat UHR_Rep2.bam > UHR_Rep2_flagstat.txt
samtools flagstat UHR_Rep3.bam > UHR_Rep3_flagstat.txt



# Run various picard CollectMetrics tools
module load picard/2.17.6
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectInsertSizeMetrics I=HBR_Rep1.bam O=HBR_Rep1_insert_size_metrics.txt H=HBR_Rep1_insert_size_metrics.pdf
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectInsertSizeMetrics I=HBR_Rep2.bam O=HBR_Rep2_insert_size_metrics.txt H=HBR_Rep2_insert_size_metrics.pdf
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectInsertSizeMetrics I=HBR_Rep3.bam O=HBR_Rep3_insert_size_metrics.txt H=HBR_Rep3_insert_size_metrics.pdf
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectInsertSizeMetrics I=UHR_Rep1.bam O=UHR_Rep1_insert_size_metrics.txt H=UHR_Rep1_insert_size_metrics.pdf
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectInsertSizeMetrics I=UHR_Rep2.bam O=UHR_Rep2_insert_size_metrics.txt H=UHR_Rep2_insert_size_metrics.pdf
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectInsertSizeMetrics I=UHR_Rep3.bam O=UHR_Rep3_insert_size_metrics.txt H=UHR_Rep3_insert_size_metrics.pdf

# Picard CollectAlignmentSummaryMetrics
cp /hpcdata/bio_data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/chr22.fa .
samtools faidx chr22.fa

java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics I=HBR_Rep1.bam O=HBR_Rep1_alignment_size_metrics.txt R=chr22.fa
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics I=HBR_Rep2.bam O=HBR_Rep2_alignment_size_metrics.txt R=chr22.fa
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics I=HBR_Rep3.bam O=HBR_Rep3_alignment_size_metrics.txt R=chr22.fa
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics I=UHR_Rep1.bam O=UHR_Rep1_alignment_size_metrics.txt R=chr22.fa
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics I=UHR_Rep2.bam O=UHR_Rep2_alignment_size_metrics.txt R=chr22.fa
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics I=UHR_Rep3.bam O=UHR_Rep3_alignment_size_metrics.txt R=chr22.fa

# multiqc
module load multiqc
multiqc .



