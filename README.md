![IGClogo](https://github.com/dsobral/ADER17S/raw/master/Logo_IGC_2014.png "IGC")
![ELIXIRPTlogo](https://github.com/dsobral/ADER17S/raw/master/elixirportugal-logo.png "ELIXIR Portugal")


# ADER17S #

# Analysis of Differential Expression with RNAseq #

A hands-on training course at Instituto Gulbenkian de Ciência (4-days)

Official course page of the Gulbenkian Training Programme in Bioinformatics - GTPB

http://gtpb.igc.gulbenkian.pt/bicourses/ADER17S/



# Overview

High-throughput technologies allow us to detect transcripts present in a cell or tissue. This introductory course covers practical aspects of the analysis of differential gene expression by RNAseq. Participants will be presented with real world examples and work with them in the training room, covering all the steps of RNAseq analysis, from planning the gathering of sequence data to the generation of tables of differentially expressed gene lists and visualization of results. We we will also cover some of the initial steps of secondary analysis, such as functional enrichment of the obtained gene lists.


# Target Audiences

Life Scientists who want to be able to use NGS data to evaluate gene expression (RNAseq). Computational researchers that wish to get acquainted with the concepts and methodologies used in RNAseq are also welcome. 


# Pre-requisites

Familiarity with elementary statistics and a few basics of scripting in R will be helpful.

Please have a look at the following resources and gauge your ability to use R in statitics at the basic level: [Coursera videos](http://blog.revolutionanalytics.com/2012/12/coursera-videos.html); [Introduction to r](http://bitesizebio.com/webinar/20600/beginners-introduction-to-r-statistical-software)

Basic Unix command line skills, such as being able to navigate in a directory tree and copy files. See, for example, ["Session 1" of the Software Carpentry training for a Unix introduction](http://bioinformatics-core-shared-training.github.io/shell-novice/). 


# Learning Objectives

Course participants will go through a series of experiences that utimately lead to create enhanced capabilities to:

1. List broad characteristics of NGS technologies and choose adequate sequencing for your biological question
2. Have a broad overview of the steps in the analysis of RNA-Seq differential expression experiments
3. Assess the general quality of the raw data from the sequencing facility
4. Do simple processing operations in the raw data to improve its quality
5. Generate alignments against a reference genome
6. Assess the general quality of the alignments and detect possible problems
7. Generate tables of counts using the alignment and a reference gene annotation
8. Generate lists of differentially expressed genes, at least for a simple pairwise comparison
9. Perform simple functional enrichment analysis and understand the concepts behind them

For this, we are providing small example datasets and exercises that participants can use to learn. 


## Learning outcomes (LO) for each unit:

### LO 1 - Plan your experiment using NGS technologies:

####	LO 1.1 - List possibilities and limitations of NGS sequencing technologies
		What choices do you have when sending your samples to the sequencing facility

####	LO 1.2 - Choose adequate sequencing for your biological question
		How do the sequencing choices influence the kind of questions you can answer

### LO 2 - List steps in the analysis of RNA-Seq differential expression experiments
		What are the steps in RNA-Seq data analysis

### LO 3 - Assess the general quality of the raw data from the sequencing facility

#### 	LO 3.1 - Interpret what are fastq files and what is their content
		What information is in fastq files, and how is it organized

#### 	LO 3.2 - Use software like FastQC to process fastq files and produce QC reports

#### 	LO 3.3 - Read QC reports of raw data to assess the general quality of data and presence of sequence bias
		Detect low quality bases in the QC reports
		Detect sequence bias and possible presence of adaptors and other contaminants

### LO 4 - Do simple processing operations in the raw data to improve its quality

#### 	LO 4.1 - Use tools such as seqtk and trimmomatic to remove low quality bases from your reads
		Use seqtk to remove a fixed number of bases from either ends of a fastq
		Use seqtk to remove low quality bases from end of a fastq file
		Use trimmomatic to filter/trim low quality bases using more complex approaches

#### 	LO 4.2 - Use tools such as cutadapt and trimmomatic to remove adaptors and other artefactual sequences from your reads
		Remove Illumina adaptor from an example dataset using cutadapt
		Remove PolyA from an example dataset using cutadapt
		Check results using FastQC on filtered data

### LO 5 - Generate alignments of processed reads against a reference genome

#### 	LO 5.1 - What is a reference genome, versioning and where to obtain genomes
		Are genomes constant?
		Obtain genome fasta from Ensembl

#### 	LO 5.2 - Alignment software: tophat2/hisat2; bwa; sailfish/salmon
		What are the conditions of using burrows-wheeler approaches?	
		Prepare a reference genome to use with hisat2 and bwa

#### 	LO 5.3 - Run an alignment: the SAM/BAM alignment format
		Run hisat2 / bwa mem in an example dataset
		What is the SAM format; what is the BAM format
	
### LO 6 - Assess the general quality of the alignments and detect possible problems

#### 	LO 6.1 - What is a reference gene annotation, versioning and where to obtain
		What is the GFF/GTF format
		Obtain genome GTF from Ensembl
			
#### 	LO 6.2 - Visualizing alignments in IGV for single genes

#### 	LO 6.3 - Use tools such as RSeQC and Qualimap to assess quality of alignments
		Interpret general alignment statistics such as percentage of aligned reads
		Check the reports to assess RNA integrity and diversity

### LO 7 - Generate tables of counts using the alignment and a reference gene annotation

#### 	LO 7.1 - The process of generating gene counts from genome aligments
		What parameters we need to consider when counting

#### 	LO 7.2 - Use tools such as htseq-counts and featurecounts to generate table of gene counts

### LO 8 - Generate lists of differentially expressed genes, at least for a simple pairwise comparison

#### 	LO 8.1 - Using the R package edgeR and DESeq2 to produce a pairwise differential expression analysis
		Use Galaxy to produce differentially expressed genes with edgeR and DESeq2
		Use edgeR and DESeq2 in R and RStudio 

#### 	LO 8.2 - Interpretation and visualization of results
		Produce PCA plots comparing all samples: outlier detection
		Visualize expression profiles of top differentially expressed genes
		Produce other plots such as vulcano plots

#### 	LO 8.3 - Use more complex settings: Generalized Linear Models
		Account for confounders using Generalized Linear Models
		Performing ANOVA-like comparisons

### LO 9 - Perform simple functional enrichment analysis and understand the concepts involved
		
#### 	LO 9.1 - How to extract meaning from a list of genes
		What are functional annotations, what types exist, and where to get them
       
#### 	LO 9.2 - Understand the concept of functional enrichment analysis, and the statistics involved
		When and why do we need multiple test corrections

#### 	LO 9.3 - Interpreting the results of functional enrichment analysis
		Using functional enrichment analysis with your lists of genes


## Detailed Program

TODO: Links to the material

### Monday, December 4th

+ 09:30 - 10:00 Introduction to the course and self presentation of the participants
+ 10:00 - 11:00 [Possibilities and limitations of NGS sequencing technologies. Choose adequate sequencing for your biological question](material/ADER17_practical.md#LO1)
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 [Steps in the analysis of RNA-Seq differential expression experiments](material/ADER17_practical.md#LO2)
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 [Interpret what are fastq files and what is their content](material/ADER17_practical.md#LO3.1). [Use software like FastQC to process fastq files and produce QC reports](material/ADER17_practical.md#LO3.2). Read QC reports of raw data to assess the general quality of data and presence of sequence bias. Use tools such as seqtk, cutadapt and trimmomatic to remove low quality bases, adaptors and other artefactual sequences from your reads.
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 What is a reference genome, versioning and where to obtain genomes. Alignment software: hisat2; bwa; salmon. Run an alignment: the SAM/BAM alignment format.

### Tuesday, December 5th

+ 09:30 - 10:00 Morning Wrap-up (what have we done so far?)
+ 10:00 - 11:00 Visualizing alignments in IGV for single genes.
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 Use tools such as RSeQC and Qualimap to assess quality of alignments
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 What is a reference gene annotation, versioning and where to obtain. The process of generating gene counts from genome aligments. Use tools such as htseq-counts to generate table of gene counts.
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 Using the R package edgeR and DESeq2 in Galaxy to produce a pairwise differential expression analysis


### Wednesday, December 6th

+ 09:30 - 10:00 Morning Wrap-up (what have we done so far?)
+ 10:00 - 11:00 Use edgeR and DESeq2 in R and RStudio. 
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 Interpretation and visualization of results.
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 Interpretation and visualization of results.
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 Use more complex settings: Generalized Linear Models


### Thursday, December 7th

+ 09:30 - 10:00 Morning wrap-up (what have we done so far?)
+ 10:00 - 11:00 Use more complex settings: Generalized Linear Models
+ 11:00 - 11:30 ***Coffee Break***
+ 11:30 - 12:30 Functional annotations: what are these and where to get them. The statistics behind functional enrichment analysis.
+ 12:30 - 14:00 ***LUNCH BREAK***
+ 14:00 - 16:00 Using functional enrichment analysis with your lists of genes.
+ 16:00 - 16:30 ***Tea Break***
+ 16:30 - 18:00 Final wrap-up Session.

