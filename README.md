# RNASeq QC and Alignment using NIAID Locus HPC
#### CREATED: 05/18/2020 
* Introductory Practical Training
* The material shown below has been adapted from open access resources listed in the References section below.   We appreciate the various groups that made high quality course material publicly available.
* The adaptation for running this data analysis using the NIAID Locus cluster has been done by **Mariam Quiñones. 

***
For an introductory lesson on using the Locus HPC - see the following link created by Dr. Poorani Subramanian [Locus](https://github.com/niaid/NGS_Intro/blob/master/notes/locus.md) 

---

## Learning Objectives:

* Describe and implement the RNA-seq workflow from the pre-processing step to the alignment of reads to the reference genome 
* Describe tools and methods within the RNA-seq workflow
* Assessing input and output filetypes
* Learn how to view alignment files on a genome browser

## Setting up to run the RNA-seq workflow

To get started with this lesson, we will start an interactive session and ask for 6 cores, by adding `pe threaded 6` to the qsub command:

```bash
$ ssh ai-submit1.niaid.nih.gov
$ qrsh -pe threaded 6 -l quick	
```

Copy the `rnaseq_lesson1` directory to your working directory:

```bash
$ cd /hpcdata/{group}/{workingdir} 
$ cp -r /hpcdata/scratch/rnaseq_lesson1 .    # if on Helix/Biowulf - use $ cp -r /scratch/rnaseq_lesson1 .
$ cd rnaseq_lesson1
```

Once you copy the directly, verify that you have a directory tree setup similar to that shown below. It is best practice to have all files you intend on using for your workflow present within the same directory.

```
rnaseq_lesson1/
	├── raw_data/
	   └── UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq
	   └── UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq
	   └── UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq
	   └── HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq
	   └── HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq
	   └── HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq
	   └── inputIDs.txt
	├── reference_data/
	   └── chr22.fa
	   └── genes.gtf
	   └── chr22index/
	├── results/
	├── scripts/
```

Below is a general overview of the steps involved in RNA-seq analysis.

<img src="https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/raw/master/img/RNAseqWorkflow.png" width="400">


So let's get started by listing the modules neeed for tools we will use to perform QC, trim, align and inspect the alignment.  We will load them as we need them, then unload

```bash
# Which modules to use?  Load and unload as needed

# module load STAR/2.7.3a-goolf-1.7.20 
# module load samtools
# module load fastqc/0.11.8-Java-1.8.0_45
# module load multiqc/1.7-Python-3.5.5
# module load cutadapt/2.3-foss-2016b-Python-3.7.3
# module load picard/2.17.6
# module load HTSeq/0.9.1-goolf-1.7.20-Python-2.7.9
```

Let's  quickly inspect the quality of input files using FastQC as learned in previous session
```bash
module load fastqc/0.11.8-Java-1.8.0_45
module load multiqc/1.7-Python-3.5.5

cd raw_data
fastqc *fastq
multiqc .

module unload fastqc/0.11.8-Java-1.8.0_45
module unload multiqc/1.7-Python-3.5.5
```

## Read Trimming
## Quality Control (*Optional*) - Trimming 

We want to make sure that as many reads as possible map or align accurately to the genome. To ensure accuracy, only a small number of mismatches between the read sequence and the genome sequence are allowed, and any read with more than a few mismatches will be marked as being unaligned. 

Therefore, to make sure that all the reads in the dataset have a chance to map/align to the genome, unwanted information can be trimmed off from every read, one read at a time. The types of unwanted information can include one or more of the following:
- leftover adapter sequences
- known contaminants (strings of As/Ts, other sequences)
- poor quality bases at read ends
- bases with bias at the 5' end of read

**We will not be performing an adapter or low quality trimming step** because:
* our data does not have an appreciable amount of leftover adapter sequences or other contaminating sequences based on FastQC.
* the alignment tool we have picked (STAR) is able to account for low-quality bases at the ends of reads when matching them to the genome.

If you need to perform trimming on your fastq data to remove unwanted sequences/bases, a recommended tool is [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html) but Btrim and Trimmomatic are also very good and fast.   A typical use of cutadapt is to remove adapter using a command line such as: $ cutadapt --adapter=AGATCGGAAGAG --minimum-length=25  -o myfile_trimmed.fastq.gz myfile.fastq.gz  

**Exercise 1**
* For this session, we will use cutadapt to remove the first 10 bases of each read.  We will use a loop to trim all files at once.  For this loop, we had created a text file with part of the IDs called "inputfastqID.txt" using the command line $ ls *fastq.gz | sed 's/\.read[1-2].fastq.gz.*//g' | sort | uniq > inputfastqID.txt 

```bash
module load cutadapt/2.3-foss-2016b-Python-3.7.3
mkdir trimmedreads

file="inputIDs.txt"
while read line
    do
    cat $line | cutadapt -u 10 -U 10 --minimum-length=25 -o trimmedreads/${line}.read1.trimmed.fastq -p trimmedreads/${line}.read2.trimmed.fastq ${line}_Build37-ErccTranscripts-chr22.read1.fastq ${line}_Build37-ErccTranscripts-chr22.read2.fastq
    echo "$line"
done <"$file"
module unload cutadapt/2.3-foss-2016b-Python-3.7.3
```

After trimming, cutadapt can remove any reads that are too short to ensure that you do not get spurious mapping of very short sequences to multiple locations on the genome. In addition to the common adapter trimming step, cutadapt can trim off any low-quality bases too, but **please note that quality-based trimming is not considered best practice, since majority of the newer, recommended alignment tools can account for this.**

## Inpection after trimming
In order to verify that cutadapt did a good job trimming, rerun **fastqc** and **multiqc**

```bash
module load fastqc/0.11.8-Java-1.8.0_45
module load multiqc/1.7-Python-3.5.5

cd trimmedreads
fastqc *trimmed.fastq
multiqc -f .

module unload fastqc/0.11.8-Java-1.8.0_45
module unload multiqc/1.7-Python-3.5.5
```

## Read Alignment
The alignment process consists of choosing an appropriate reference genome to map our reads against, and performing the read alignment using one of several splice-aware alignment tools such as [STAR](https://github.com/alexdobin/STAR) or [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) (HISAT2 is a successor to both HISAT and TopHat2). The choice of aligner is a personal preference and also dependent on the computational resources that are available to you.
 
For this workshop we will be using STAR (Spliced Transcripts Alignment to a Reference), an aligner designed to specifically address many of the challenges of RNAseq read mapping. 

### STAR Alignment Strategy

STAR is shown to have **high accuracy** and outperforms other aligners by more than a **factor of 50 in mapping speed (but also requires quite a bit of memory)**. The algorithm achieves this highly efficient mapping by performing a two-step process:

1. Seed searching
2. Clustering, stitching, and scoring

#### Seed searching

For every read that STAR aligns, STAR will search for the longest sequence that exactly matches the reference genome:

![STAR_step1](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/img/alignment_STAR_step1.png)
	
The different parts of the read that are mapped separately are called 'seeds'. So the first MMP that is mapped to the genome is called *seed1*.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, which will be *seed2*. 

![STAR_step2](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/img/alignment_STAR_step2.png)

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the longest matching portions of the read, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping. More details on the algorithm itself can be found in the [STAR publication](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635). 

**If STAR does not find an exact matching sequence** for each part of the read due to mismatches or indels, the seed will be extended.

![STAR_step3](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/img/alignment_STAR_step3.png)

**If extension does not give a good alignment**, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

![STAR_step4](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/img/alignment_STAR_step4.png)


#### Clustering, stitching, and scoring

The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of seeds that have good alignment scores and are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.). 

![STAR_step5](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/img/alignment_STAR_step5.png)

### Running STAR

Aligning reads using STAR is a two-step process:   

1. Create a genome index 
2. Map reads to the genome

> A quick note on shared databases for human and other commonly used model organisms. The Locus cluster has a designated directory at `/hpcdata/bio_data/` in which there are files that can be accessed by any user. These files contain, but are not limited to, genome indices for various tools, reference sequences, tool specific data, and data from public databases, such as NCBI. So when using a tool and requires a reference of sorts, it is worth taking a quick look here because chances are it's already been taken care of for you. 

> If the desired reference is not already available on the bio_data directory, contact the locus team and request it.  Also, reference genomes can be downloaded from various sources (e.g. http://hgdownload.soe.ucsc.edu/downloads.html#human)

```bash
$ ls -l /hpcdata/bio_data/iGenomes/
```

#### Creating a genome index

* For this training we are using reads that originate from chr22 therefore we will only need a genome index for chr22.  Since it is a small chromosome, the step to generate the index should only take a couple of minutes. 

* Please remember that for a regular analysis, you should make sure to use an existing index for the entire genome or create your own, which will take a while.  For indexing the reference genome, a reference genome (FASTA) is required and an annotation file (GTF or GFF3) is suggested for a more accurate alignment of the reads. For our exercise, the FASTA was copied from /hpcdata/bio_data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/chr22.fa) and the annotations file copied from /hpcdata/bio_data/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf

The basic options to **generate genome indices** using STAR are as follows:

* `--runThreadN`: number of threads
* `--runMode`: genomeGenerate mode
* `--genomeDir`: /path/to/store/genome_indices
* `--genomeFastaFiles`: /path/to/FASTA_file (reference genome)
* `--sjdbGTFfile`: /path/to/GTF_file (gene annotation)
* `--sjdbOverhang`: readlength -1

```bash
module load STAR/2.7.3a-goolf-1.7.20 

cd ../../ # you should now be in the directory `rnaseq_lesson1`
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir reference_data/chr22index \
--genomeFastaFiles reference_data/chr22.fa \
--sjdbGTFfile reference_data/genes.gtf \
--sjdbOverhang 89 \
--genomeSAindexNbases 11
```

The basic options for **mapping reads** to the genome using STAR are as follows:

* `--runThreadN`: number of threads
* `--readFilesIn`: /path/to/FASTQ_file
* `--genomeDir`: /path/to/genome_indices
* `--outFileNamePrefix`: prefix for all output files

We will also be using some advanced options:
* `--outSAMtype`: output filetype (SAM default)
* `--outSAMUnmapped`: what to do with unmapped reads
* `--outSAMattributes`: SAM attributes

Note that default filtering is applied in which the maximum number of multiple alignments allowed for a read is set to 10. If a read exceeds this number there is no alignment output. To change the default you can use `--outFilterMultimapNmax`, but for this lesson we will leave it as default. The advanced parameters that we are going to use are described below:

More details on STAR and its functionality can be found in the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), we encourage you to peruse through to get familiar with all available options.

Now let's put it all together! The full STAR alignment command is provided below.

> If you like you can copy-paste it directly into your terminal. Alternatively, you can manually enter the command, but it is advisable to first type out the full command in a text editor (i.e. [Sublime Text](http://www.sublimetext.com/) or [Notepad++](https://notepad-plus-plus.org/)) on your local machine and then copy paste into the terminal. This will make it easier to catch typos and make appropriate changes. 

## Let's align one pair of trimmed reads to the indexed genome (chr22)
* Where do the input files come from?  
* See additional details: https://github.com/griffithlab/rnaseq_tutorial/wiki/RNAseq-Data .  Paired-end 101-mers generated on an Illumina HiSeq instrument.  UHR (Universal Human Reference), HBR (Human Brain Reference).  These datasets are highlighted on the tutorial https://github.com/griffithlab/rnaseq_tutorial/wiki/Differential-Expression
* UHR + ERCC Spike-In Mix1, Replicate 1
* UHR + ERCC Spike-In Mix1, Replicate 2
* UHR + ERCC Spike-In Mix1, Replicate 3
* HBR + ERCC Spike-In Mix2, Replicate 1
* HBR + ERCC Spike-In Mix2, Replicate 2
* HBR + ERCC Spike-In Mix2, Replicate 3

**Exercise 2**
* Align one pair to chr22 index
Make sure you are in /rnaseq_lesson1/raw_data/
```bash
# do HBR_Rep1
STAR --genomeDir ../reference_data/chr22index \
--runThreadN 6 \
--readFilesCommand trimmedreads/HBR_Rep1_ERCC-Mix2.read1.trimmed.fastq HBR_Rep1_ERCC-Mix2.read2.trimmed.fastq \
--outFileNamePrefix results/HBR_Rep1_test \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 

# or UHR_Rep3
STAR --genomeDir ../reference_data/chr22index \
--runThreadN 6 \
--readFilesCommand trimmedreads/UHR_Rep3_ERCC-Mix1.read1.trimmed.fastq trimmedreads/UHR_Rep3_ERCC-Mix1.read2.trimmed.fastq \
--outFileNamePrefix results/UHR_Rep3_test \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard
```

* Using the `less` command take a look at `*.final.out` and answer the following questions:  
	1. How many reads are uniquely mapped?
	2. How many reads map to more than 10 locations on the genome?
	3. How many reads are unmapped due to read length?
	
* Also, if desired use samtools flagstat to view percentage of mapped reads and other statistics
```bash
module load samtools
samtools flagstat name_of_bamfile.bam
modulle unload samtools
```


### Please know that on a real analysis, you should align all samples.  If needed, merge technical replicates as appropriate and then visualize on a browser.  For alignment of many samples, use a loop as shown below for convenience.
***DO NOT WORRY ABOUT RUNNING THE CODE BELOW DURING THIS LIVE TRAINING SESSION*** Instead after this class, try running the entire script that is available on the scripts directory.  For the purpose of this training, we will provide the results below for you to use in subsequent steps. 

```bash
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
```

### Alignment Outputs (SAM/BAM)

The output we requested from STAR is a BAM file, and by default returns a file in SAM format. **BAM is a binary version of the SAM file, also known as Sequence Alignment Map format.** The SAM file is a tab-delimited text file that contains information for each individual read and its alignment to the genome. The file begins with an optional header (which starts with '@'), followed by an alignment section in which each line corresponds to alignment information for a single read. **Each alignment line has 11 mandatory fields** for essential mapping information and a variable number of fields for aligner specific information.

These fields are described briefly below, but for more detailed information the paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) is a good start.

![SAM](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/raw/master/img/sam_bam.png)

![SAM](https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/raw/master/img/sam_bam3.png)


* The SAM flags are explained in more details here: [picard_explain-flags](https://broadinstitute.github.io/picard/explain-flags.html)

* We are now getting ready to view the alignment files (BAM) in the IGV browser.  Before doing so, we should rename files, merge and index

* If you were not able to align the reads during the class, jump to the **Merge alignment** section where you will find a copy of the aligned bam files.  On the other hand, if you have bam files resulting from the alignment exercise using STAR, proceed to rename the files. 

```bash
cd ../results
cp HBR_Rep1*bam HBR_Rep1.bam
cp HBR_Rep2*bam HBR_Rep2.bam
cp HBR_Rep3*bam HBR_Rep3.bam
cp UHR_Rep1*bam UHR_Rep1.bam
cp UHR_Rep2*bam UHR_Rep2.bam
cp UHR_Rep3*bam UHR_Rep3.bam
```

```bash
# optional inspection
samtools view -h results/HBR_Rep1.bam | less
```
Scroll through the SAM file and see how the fields correspond to what we expected.

View statistics such as % aligned using samtools flagstat
```bash
samtools flagstat results/HBR_Rep1.bam
```
### Merge alignment
**Exercise 3: Merge, Index and View files using IGV browser**

```bash

# If you have not yet aligned all 6 samples, copy the directory rnaseq_bams_lesson1 to the working directory
cp -r /hpcdata/scratch/rnaseq_bams_lesson1 .    # if on Biowulf/Helix, use cp -r /scratch/rnaseq_bams_lesson1 .

module load picard/2.17.6

cd rnaseq_bams_lesson1
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar MergeSamFiles OUTPUT=UHR.bam INPUT=UHR_Rep1.bam INPUT=UHR_Rep2.bam INPUT=UHR_Rep3.bam
java -Xmx2g -jar ${EBROOTPICARD}/picard.jar MergeSamFiles OUTPUT=HBR.bam INPUT=HBR_Rep1.bam INPUT=HBR_Rep2.bam INPUT=HBR_Rep3.bam

module unload picard/2.17.6
```

### Index bam files
```bash
module load samtools

# create the index
samtools index UHR.bam
samtools index HBR.bam

module unload samtools
```
### Load bam files to the IGV browswer

Use _**locus mounted drive** to upload files to IGV:
smb://locusfileserver.niaid.nih.gov/{group}/{directory}/

> **NOTE: You can also transfer files to your laptop using the command line**
>
> Similar to the `cp` command, there is a command that allows you to securely copy files between computers. The command is called `scp` and allows files to be copied to, from, or between different hosts. It uses ssh for data transfer and provides the same authentication and same level of security as `ssh`. 
>
> First, identify the location of the _origin file_ you intend to copy, followed by the _destination_ of that file. Since the original file is located on Orchestra, this requires you to provide remote host and login information.

> ```bash
> $ scp user_name@ai-submit1.niaid.nih.gov:/hpcdata/group/directory/rnaseq_lesson1/results/*.bam /path/to/directory_on_laptop
> ```

**Visualize**

* Start [IGV](https://www.broadinstitute.org/software/igv/download), *you should have this previously installed on your laptop*.
* Load the Human genome (hg19) into IGV using the dropdown menu at the top left of your screen. 
**Note**: there is also an option to "Load Genomes from File..." under the "Genomes" pull-down menu - this is useful when working with non-model organisms.  Select chr22.
* Load the .bam file using the **"Load from File..."** option under the **"File"** pull-down menu. *IGV requires the `.bai` file to be in the same location as corresponding `.bam` file that you want to load into IGV, but there is no other direct use for this index file.*

* Explore the read alignment in the following genes: EIF3L, NDUFA6, and RBX1
* Also, noticed that SULT4A1 and GTSE1 are differentially expressed. Are they up-regulated or down-regulated in the brain (HBR) compared to cancer cell lines (UHR)?

* The browser should look like this image
![GTSE1](https://github.com/niaid/RNASeqQC_Alignment/blob/master/GTSE1.png)

***

#### ***Challenge - Post Alignment QC***
> For visualizing metrics such as insert size, run Picard tools and FastQC using the script [alignment_postQC.sh](https://github.com/niaid/RNASeqQC_Alignment/blob/master/alignment_postQC.sh) which is a short version of the example here:
https://pmbio.org/module-03-align/0003/05/01/PostAlignment_QC/

***
#### Optional steps for obataining counts from the alignment files
**Exercise 4: use HTSeq for generating count table
* Run htseq-count on alignments instead to produce raw counts instead of FPKM/TPM values for differential expression analysis. Refer to the HTSeq documentation for a more detailed explanation:http://www-huber.embl.de/users/anders/HTSeq/doc/count.html

* Below are the steps for generating a table with counts for each gene per sample.  After the count table, feel free to try a differential expression workflow as described [here](https://github.com/griffithlab/rnaseq_tutorial/wiki/Differential-Expression). 
```bash
module load HTSeq/0.9.1-goolf-1.7.20-Python-2.7.9


htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id UHR_Rep1.bam ../reference_data/genes.gtf > UHR_Rep1_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id UHR_Rep2.bam ../reference_data/genes.gtf > UHR_Rep2_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id UHR_Rep3.bam ../reference_data/genes.gtf > UHR_Rep3_gene.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id HBR_Rep1.bam ../reference_data/genes.gtf > HBR_Rep1_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id HBR_Rep2.bam ../reference_data/genes.gtf > HBR_Rep2_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id HBR_Rep3.bam ../reference_data/genes.gtf > HBR_Rep3_gene.tsv
```

* Merge results files into a single matrix for use in edgeR. The following joins the results for each replicate together, adds a header, reformats the result as a tab delimited file, and shows you the first 10 lines of the resulting file :

```bash
join UHR_Rep1_gene.tsv UHR_Rep2_gene.tsv | join - UHR_Rep3_gene.tsv | join - HBR_Rep1_gene.tsv | join - HBR_Rep2_gene.tsv | join - HBR_Rep3_gene.tsv > gene_read_counts_table_all.tsv
echo "GeneID UHR_Rep1 UHR_Rep2 UHR_Rep3 HBR_Rep1 HBR_Rep2 HBR_Rep3" > header.txt
cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > gene_read_counts_table_all_final.tsv
rm -f gene_read_counts_table_all.tsv header.txt
head gene_read_counts_table_all_final.tsv
```


## References
1. "RNA-Seq workflow" (https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/lessons/07_rnaseq_workflow.md) author: "Mary Piper, Meeta Mistry, Radhika Khetani,  Bob Freeman". date: "Tuesday, August 22, 2017". [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/).  Alignment step - https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html and Introductory slides - https://github.com/hbctraining/Intro-to-rnaseq-hpc-orchestra/blob/master/lectures/rna-seq_design.pdf
2. Description on QC errors - https://www.epigenesys.eu/images/stories/protocols/pdf/20150303161357_p67.pdf
3. STAR publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/ and Manual https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
4. RNAseq tutorial: https://github.com/griffithlab/rnaseq_tutorial/wiki/RNAseq-Data and https://rnabio.org/
5. Malachi Griffith*, Jason R. Walker, Nicholas C. Spies, Benjamin J. Ainscough, Obi L. Griffith*. 2015. Informatics for RNA-seq: A web resource for analysis on the cloud. PLoS Comp Biol. 11(8):e1004393
6. QC after alignment - http://rseqc.sourceforge.net/ and example: https://github.com/griffithlab/rnaseq_tutorial/wiki/PostAlignment-QC
7. HTSEq https://htseq.readthedocs.io/en/release_0.9.1/count.html


***
*This lesson shared at https://github.com/hbctraining/ was developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
*In addition, the materials used in this lesson were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
