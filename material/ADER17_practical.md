### Learning Outcome 1: Plan your experiment using NGS technologies

A good source of information for this part is [RNA-seqlopedia](http://rnaseq.uoregon.edu).

#### Obtaining your RNA

The first step in a transcriptomic experiment is to obtain the RNA. After isolating total RNA from cells, one can directly sequence it. Nonetheless, the majority of the RNA in a cell is ribosomal RNA, which may need to be removed using specific kits. Moreover, total RNA also contains unprocessed immature transcripts and RNA targeted for degradation (at different stages of processing). 

Therefore, unless one is interested in non-coding RNAs or other aspects related to transcription, it is usually better to apply protocols that extract the mature mRNAs (usually through the PolyA tails). Since most people are interested in coding-genes, it is more common to use mRNA-specific protocols. 

Some protocols can also keep strand information. In this case, the reads have the same (or the reverse) strand as the transcribed RNA. This is particularly relevant when sequencing total RNA, noticeably to distinguish real transcripts from transcriptional activity resulting from stalled promoters or enhancers. It can also be useful to distinguish between overlapping genes.
 
Finally, we also need to consider the amount of material available. Are we dealing with samples with a lot of RNA (eg. cell cultures), or short amounts (eg. small tissue samples, single-cell) that are prone to amplification artifacts and presence of contaminant sequences? 

All these aspects need to be taken into consideration when analyzing the data.

#### <a id="LO1">Options for sequencing</a>

At the moment, the sequencing technology most often used (by far) is Illumina. The following links are a good source of information regarding this sequencing technology:
* [Illumina Sequencing by Synthesis](https://www.youtube.com/watch?&v=fCd6B5HRaZ8).
* [Elaine Mardis talk on NGS](https://www.youtube.com/watch?v=v1DbcJD4Ry0).

Options to be considered when sequencing:
* Single versus Paired-end
* Read Length
* Coverage (number of reads)

For the analysis of differential gene expression, long reads, paired-end, and stranded library preparation methods are not as important, particularly if a reference genome is available. Focus should be given on replicates in order to obtain accurate measures of variances. The number of replicates and depth of sequencing depends on the experiment. For highly controlled conditions (such as cell cultures), 2-3 replicates could be enough. In terms of coverage, 10-40M reads should be enough to capture most "reasonably" expressed genes, although in single-cell experiments less reads may be sufficient (5-10M). Nonetheless, to be able to more accurately estimate how much is needed, one should always generate small pilot datasets. 

For this course, we will focus on the analysis of differential gene expression between two conditions. Thus, we assume unstranded mRNA-specific library preparation methods, sequenced using illumina (NextSeq, HiSeq) short (less than 100bp) single-end reads. We also assume 2-3 replicates per condition, sequenced to a medium throughput (10-40M reads). We will nonetheless briefly discuss what to do in other cases such as longer reads, paired data, stranded data, and more complex differential expression conditions. 

### Learning Outcome 2: List steps in the analysis of RNA-Seq differential expression experiments

Steps in the analysis of RNA-Seq:
* QC of Raw Data; (Learning Outcome 3)
* Preprocessing of Raw Data (if needed); (Learning Outcome 4) 
* Alignment of “clean” reads to reference genome (Learning Outcome 5)
* QC of Aligments (Learning Outcome 6)
* Generate table of counts of genes/transcripts (Learning Outcome 7)
* Differential Analysis tests (Learning Outcome 8)
* Post-analysis: Functional Enrichment (Learning Outcome 9)

# Learning Outcome 3: Assess the general quality of the raw data from the sequencing facility

## LO 3.1 - Interpret what are fastq files and what is their content

Most high-throughput sequencing (HTS) machines output [fastq files](https://en.wikipedia.org/wiki/FASTQ_format), the “de facto” current standard in HTS. Fastq files are simply text files, where each block of information (a sequenced DNA fragment, or read) in this format is encoded as 4 lines:

	@read_identifier
	read_sequence
	+ separator line
	base_qualities

Each base has a quality character associated with it, representing how confidently the machine identified (called) the base. The probability of error per base is given as a [Phred score](https://en.wikipedia.org/wiki/Phred_quality_score), calculated from an integer value (Q) derived from the quality character associated to the base. Useful reference values of Q include:
* Q=10 - 90% accuracy
* Q=20 - 99% accuracy
* Q=30 - 99.9% accuracy
* Q=40 - 99.99% accuracy

Although there's theoretically no limit, Q usually goes up to around 40 in recent illumina machines.

You can see a few fastq files in the folder fastq_examples:
* sample_quality_and_adaptors.fastq.gz
* sample_adaptors.fastq.gz
* 20150821.A-2_BGVR_P218_R1.sample.fastq.gz
* 20150821.A-2_BGVR_P218_R2.sample.fastq.gz

Since each fastq can have several million reads, they can become very big. Therefore, it is usual to keep them in a compressed format such as gzip. Most recent software can directly read compressed fastq files.

**Task**: Upload all sample files into your Galaxy and inspect them

You probably noticed that two of the example files have the same name, except for R1 and R2. This is an example of a paired-end dataset. If you inspect both datasets, you can find the same identifiers in each of the files, in the same order. In R1 you have the forward reading of a fragment, and in R2 you have the reverse reading of the same fragment.


## LO 3.2 - Use software like FastQC to process fastq files and produce QC reports

High Throughput Sequencing machines read thousands or millions of sequences in paralell. As you can imagine, this usually generates large fastq files, with millions of lines. Manually inspecting quality of each read is out of question. Specialized software has been developed to provide quality measures for fastq files generated by HTS machines. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular program to generate quality reports on fastq data. Running FastQC on your raw data is usually the first thing you should do once you receive a new dataset.

**Task**: In Galaxy, run FastQC in each of the example files.


## LO 3.3 - Read QC reports of raw data to assess the general quality of data and presence of sequence bias

FastQC reports provide a series of plots that allow the user to assess the overall quality of their raw data and detect potential biases and problems. 

One of the plots indicates distribution of base qualities along the length of reads. You can notice that, at least for illumina data, on average the quality of each base tends to decrease along the length of the read. You can also see that the reverse read (R2) is usually of worse quality than the forward read (R1). Therefore, short single-end reads usually have better average quality, and are often ready to use right out of the sequencer.

![Base Quality](images/base_quality.png) ![Tile Quality](images/tile_quality.png)

Other plots indicate biases in nucleotidic content of reads, either globally (as %GC plots), or positionally. Global bias in nucleotidic content can be useful to search for signs of contaminants. On the other hand, positional bias are useful to detect presence of artefactual sequences in your reads such as adaptors. Another insight you may obtain from this information are potential biases in the preparation of your library. For example, random hexamer priming is actually not truly random, and preferentially selects certain sequences. The currently popular transposase-based enzymatic protocol, although reasonably random, is also not completely random, and you can see this through positional bias, particularly in the beginning of reads. The presence of adaptors is a relatively common event, and therefore specific plots exist to detect the presence of the most commonly used adaptors. Finally, the presence of repetitive sequences can also suggest contaminants, pcr artifacts, or other types of bias.

![Base Bias](images/base_bias.png) ![Adaptor](images/adaptor.png)

**Task**: Inspect the FastQC Reports generated previously and detect potential issues.

# Learning Outcome 4: Do simple processing operations in the raw data to improve its quality

In most cases, particularly if you're sequencing short, single-end reads, the quality of your raw data is good enough to continue without any preprocessing. In fact, if you send your sequencing to an external facility, they often do these verifications and filtering for you, and you have “clean” sequences in the end. Nonetheless, it is always better to check before proceeding. 

Sometimes things can go wrong, and you may need to do something about it. Some types of problems, like presence of contaminants, or some instances of positional bias will require to go back and redo the experiments. Other issues can be minimized. 


## LO 4.1 - Use tools such as seqtk and trimmomatic to remove low quality bases from your reads

As you may have noticed before, reads tend to lose quality towards their end, where there is a higher probability of erroneous bases being called. To avoid problems in subsequent analysis, you should remove regions of poor quality in your read, usually by trimming them from the end of reads using tools such as [seqtk](https://github.com/lh3/seqtk). 

**Question**: Even if all bases that your machine reads have a Q=20 (1% error rate), what is the probability that one 100bp read is completely correct? To answer this, consider also that all bases are read independently.

**Task**: In Galaxy, use seqtk_trimfq with different error thresholds in the example datasets. Use FastQC to evaluate the impact of the procedure. Compare this with the simpler approach of cutting your reads to a fixed length.

**Question**: If you are too stringent, you may remove too many bases, but if you are too lenient, you may fall in local optima, because behind a good quality base may be more bad quality ones. What other strategies you can imagine to filter your reads?

Another popular tool to filter fastq files is [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). This tool implements more ellaborate trimming strategies, such as average window threshold.

**Task**: In Galaxy, use Trimmomatic to remove low quality bases from the example datasets. Notice that the default method in Trimmomatic is a 4bp window average, with a threshold of Q=20. Finally, look at the impact using FastQC of trimmed reads. NOTE: Trimmomatic requires you to specify that you use the "standard" Phred Q scale (fastqsanger), which was different from the one used in older datasets (before 2012), so you need to manually change the datatype of your dataset from generic fastq to fastqsanger.

## LO 4.2 - Use tools such as cutadapt to remove adaptors and other artefactual sequences from your reads

Sequence machines often require that you add specific sequences (adaptors) to your DNA so that it can be sequenced. For many different reasons, such sequences may end up in your read, and you usually want to remove these adaptors. Moreover, cDNAs may contain parts of the non-genomic polyA tails that are part of mature mRNAs. Since these sequences are not part of the genome, they may prevent proper alignment and need to be removed before proceeding.

To remove these unwanted sequences, not only you have to look for the sequence in the reads, but also allow for sequencing errors, as well as the presence of incomplete sequences. Tools such as [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html) do precisely this.

**Task**: In Galaxy, use cutadapt to remove adaptors from sample_adaptors.fastq. In this sample, we know that we used the illumina adaptor GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT, so try to remove this from the 3' end of reads and see the impact of the procedure. What happened? You noticed that almost no read was affected. This is because what you get is a readthrough, so you actually have the reverse complement of the adaptor. Now, try the same procedure but with AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC (reverse complement of the previous). Much better, no? 

One issue of removing the adaptors is that you need to know which ones were used in your data. FastQC can already tell you which one was used, and you can then go to the illumina manual to search for its sequence. Since Illumina is used most of the time, these adaptors are already integrated in tools like Trimmomatic, which also take in consideration issues like reverse complement. 

**Task**: Use Trimmomatic to remove adaptors from sample_adaptors.fastq using Truseq adaptors and use FastQC to see the results.

Overall, you can use Trimmommatic to do both quality and adaptor trimming. Moreover, Trimmomatic includes widely used adaptors and also transparently handles the issue of paired-end, in case you have this. So it's probably a good choice for general use. There are, nonetheless, several alternative tools that do many of these procedures.

Task: Use Trimmommatic to do quality filtering and adaptor trimming in sample_quality_and_adaptors.fastq (use Nextera adaptors), as well as in the paired-end example RNA-Seq data (use Truseq adaptors). Use FastQC to evaluate the impact of the procedure. Notice that, if you're too strict, you may end up loosing valuable data.

**Task**: Finally, inspect a complete dataset (your own, or some we provided). First, use FastQC to detect potential issues. If necessary, use Trimmomatic (or any of the other tools we tried).

# Learning Outcome 5: Generate alignments of processed reads against a reference genome

We've checked the quality of our raw data, and did any necessary preprocessing, so we should now be ready to use it. 

## LO 5.1 - What is a reference genome, versioning and where to obtain genomes

We now need to align the reads against a reference genome. Genomes were (and are still) usually obtained through the efforts of large consortia, which eventually create portals that make the data available for the scientific community. [ENSEMBL](http://www.ensembl.org) (in Europe) and [UCSC genome browser](http://genome.ucsc.edu/) (in the US) emerged first as resources to display and explore the human data, and latter agglomerated data for other model and non-model organisms, making them very convenient resources for high quality genomes. 

Genome assemblies are continuously updated with new information, particularly for large eukaryotic genomes. Even the human genome, that was "completed" in 2001, is regularly being updated. More recent updates of the Human genome do not change the core sequence, but add for example alternative haplotypes for complex and highly variable regions such as the HLA. It is also very frequent to have several alternative genomes for the same species (eg. different lab strains of mice, or other model organisms). Moreover, large genomes contain many repetitive elements, which are usually masked for secondary analysis like gene annotation. For the alignment of NGS data, it is usually recommended to use full, unmasked, sequences. It is also common to ignore alternative haplotypes, although in human this depends on the goals of the study.

It is fundamental to register the version of the genome used, as well as from where (and when) it was obtained. When performing analysis using resources like Galaxy, genomes are often already integrated in those resources. You should always note as much information as possible about the genome you're using and, if in doubt, contact the service providers to find out more information.

Finally, another alternative is to use cDNA sequences directly as a reference. This is sometimes the only alternative, when full good quality genomes are not available. The presence of multiple alternative transcripts can make the alignment more difficult, but more recent approaches can actually take this information in consideration. We can also select collections of cDNAs that are relevant for our analysis (eg. focusing on protein-coding cDNAs, and/or choosing a single representative cDNA per gene).

**Task**: Obtain genomic fasta for Drosophila melanogaster from the Ensembl website. Finally, also download a fasta with cDNA. Take note of the Ensembl version, as well as the version of your genome (in case later you wano to integrate data that is not from Ensembl). Obtain genomic and cDNA fasta from ENSEMBL for the species relevant for your complete dataset.

## LO 5.2 - Alignment software: hisat; bwa; salmon

To be able to align millions of short reads to a (sometimes large) reference genome, novel, more efficient, alignment methods had to be developed. The most popular are based on the [burrows-wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform), of which [bwa](http://bio-bwa.sourceforge.net/) and [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) are examples. They enable alignment of millions of reads in a few minutes, even in a common laptop.  

Nonetheless, these methods rely on building large data structures that must stay entirely in memory, and thus aligning against larger genomes such as Human may require more computational resources than available in more modest laptops. Moreover, most software implementations of these methods, although open-source and freely available, are usually only available for use within a linux-based environment. In any case, given the large size of the raw data and the computational requirements for the alignment, it is most likely that at least the alignment step will have to be run in a dedicated computational server. Services such as Galaxy or Chipster hide the computational infrastructure from you and provide simple to use interfaces that allow you to run analysis that would otherwise be hard to perform on your own computer.

Methods based on the burrows-wheeler transform make assumptions to speed up the alignment process. Namely, they require the reference genome to be very similar to your sequenced DNA (less than 2-5% differences). For example, mouse data will align poorly to the human genome, although in the case of RNA-Seq this is less problematic since genes tend to be much better conserved than the rest of the genome (you would probably still bias your results to better conserved genes). Moreover, these fast alignment algorithms are not optimal, and therefore sometimes make some mistakes, although they work quite well most of the time. 

Eukaryotes contain the extra complication of splicing, where your read will be spread through multiple regions of the genome (usually, different exons of the same transcript). When using small, single-end reads, this is less of a problem, since it is less likely that your reads will overlap significantly with a splice site. Nonetheless, it is a good idea to use an aligner that allows split reads. [Hisat](https://ccb.jhu.edu/software/hisat2/index.shtml) (based on bowtie) is one such splice-aware aligner (it is an update of the better known Tophat aligner). It still uses the same approach as before, but with extensions to allow splitting of reads. Recent updates of bwa (bwa mem) also allows splitting of reads, and can be used for RNA-Seq data .

Finally, another set of more recent approaches quickly gaining in popularity align directly against the transcriptome, without the need for a reference genome. [Salmon](https://combine-lab.github.io/salmon/) provides transcript-level estimates of gene expression. These methods are very fast (mostly because they only align against the transcriptome), and use more elaborate statistical methods to handle the presence of different alternative splice forms that difficult the attribution of a read to a transcript. Some of these methods, such as salmon, also take explicitly in consideration bias related to differences in transcript length and nucleotide composition. 

## LO 5.3 - Run an alignment: the SAM/BAM alignment format

As we mentioned before, aligners for NGS data depend on large data structures for their efficiency. These structures (like the blast databases) are built from the fasta file containing the sequence of the reference genome. This process is relatively slow and computationally intensive, although it is only necessary to do it once for every reference genome. Therefore, before aligning your reads, it is necessary to do an indexing step on the genome sequence that will be used for alignment. If using the tools on the command line, one needs to explicitly perform this step. Using services such as Galaxy, this step is hidden from the user. 

When performing the alignment in Galaxy, you usually have two options: either you pass the tool a fasta with the reference genome, or you select an available genome. When using an available genome, the indexing step was already performed, while if you provide your own fasta of the genome, an indexing step will have to be performed before the alignment step. If your genome of interest is relatively large (roughly >100Mb), it is more efficient to have it pre-built, particularly if you're reusing it often. For this, you will need to ask the persons managing the service you're using.

**Task**: In Galaxy, run Hisat2 on one of the paired-end example files (in single-end mode) against the Drosophila genome that should be prebuilt in your Galaxy instance. Now run the same, but using as input the fasta for the Drosophila genome that you downloaded previously instead of the prebuilt one already available. Compare the differences in the time it takes. 

Researchers felt the need to develop new, more practical formats to store millions of alignments being generated by these aligners. The [Sequence Alignment/Map (SAM) format](https://samtools.github.io/hts-specs/SAMv1.pdf) is a tabular text file format, where each line contains information for one alignment. SAM files are most often compressed as BAM (Binary SAM) files, to reduce space and allow direct access to alignments in any arbitrary region of the genome. Several tools only work with BAM files. Some aligners still produce only SAM files, which may need to be converted to BAM.

**Task**: Upload the example SAM files in the same folder as the example fastq. Inspect them, comparing single-end to paired-end. Convert to BAM (you'll need the Drosophila genome you've uploaded).

Most genomes (particularly mamallian genomes) contain areas of low complexity, composed mostly of repetitive sequences. In the case of short reads, sometimes these align to multiple regions in the genome equally well, making it impossible to know where the fragment came from. Longer reads are needed to overcome these difficulties, or in the absence of these, paired-end data can also be used. Some aligners (such as hisat or bwa) can use information from paired reads to help disambiguate some alignments. Information on paired reads is also added to the SAM/BAM file when proper aligners are used.

**Task**: In Galaxy, run hisat2 and bwa mem with the example paired-end data. In the guilgur folder, you'll have data from [Guilgur et al, 2014](https://elifesciences.org/content/3/e02181). Run bwa mem of those files against the prebuilt Drosophila genome. 

**Task**: In Galaxy, run Hisat2 on the complete dataset. Use a genome that is already prebuilt, to save time. During the day, so that you can work on other things in Galaxy during the course, run one alignment at a time. If you still need to run alignments at the end of the day, then you can launch all in a queue to be run overnight.

Salmon directly estimates transcript expression (not alignments), and thus we will come back to it later on.

# Learning Outcome 6: Assess the general quality of the alignments and detect possible problems

## LO 6.1 - What is a reference gene annotation, versioning and where to obtain

To estimate gene expression, we need to define the genes by identifying their position and structure in the genome. This information is stored 
in a hierarchical fashion (the genes, their transcripts, each transcript's exons, and so on...) in formats such as the [Generic Feature Format (GFF) files](http://gmod.org/wiki/GFF3). These consist basically of tabular text files with positions of genes (and their components) in the genome (for a specific genome version), as well as other information about the gene such as its name. Another common format used for annotations is the [BED format](http://genome.ucsc.edu/FAQ/FAQformat.html#format1). 

Each gene annotation is deeply associated to one given version of the genome (because it contains positions in that genome), but the same genome version may (and usually has) several gene annotation versions. The same way one should keep in detail the version of the genome, we should also take note of the version of the gene annotation being used, and from where and when we obtained it.

Gene annotations are usually complex to create, particularly for large mammalian genomes, and are a permanent work in progress (even more than the genome). Annotation of the Human genes is the work of several large groups of researchers. Other model organisms (such as the mouse) also have dedicated teams to curate their genes.  Non-model organisms that are less intensively studied may suffer from having less well characterized annotations, and frequently derived from other better studied organisms. 

The same way ENSEMBL is a good source for the genome sequence, it is also a good source to obtain gene annotations. ENSEMBL even defined a specific variant of the GFF format ([GTF](http://www.ensembl.org/info/website/upload/gff.html)) which is commonly accepted by most applications. 
 
**Task**: Obtain the latest Drosophila melanogaster GTF from Ensembl, as well as the GTF for the organism relevant for your complete dataset.


## LO 6.2 - Visualizing alignments in IGV for single genes

To visualize the alignments along the reference genome one can use software such as [IGV](http://software.broadinstitute.org/software/igv/) or [Tablet](https://ics.hutton.ac.uk/tablet/), which work with the most common operating systems. To avoid loading all alignments simultaneously in memory, and to be able to quickly search for region-specific alignments, this software uses the BAM format. 

**Task**: Run IGV and look at the provided sample BAM files with alignments. In IGV, load the Drosophila genome as reference (fasta), and then load the provided annotation file (gtf) and alignment files (*.bam). Note that for this case, you need to load Drosophila BDGP5 and not BDGP6 (unless you redid alignments using the raw data that are also provided). 
- Look at position: 3L:15033260-15038204
- Look at position: X:20564838-20570348
- Look at position X:5793758-5799858

**Question**: Would you be able to detect all of what you saw here using microarrays? If not, what and why?

**Task**: Download the BAM files you generated for your complete dataset, and load in IGV. Don't forget to also download the companion bai index files. Also, don't forget you first need to load an appropriate genome of reference and gene annotation (GTF file) that you should have downloaded previously.


## LO 6.3 - Use tools such as RSeQC and Qualimap to assess quality of alignments

After generating alignments and obtaining a SAM/BAM file, how do I know this step went well? In fact, there are potential issues that we can only detect after we try to align against the reference genome. The same way FastQC generates reports of fastq files to assess quality of raw data, there are programs that generate global reports on the quality of alignments. One popular tool for this is [qualimap](http://qualimap.bioinfo.cipf.es/). 

One important general measure is how many (out of all reads) were properly aligned against the reference genome. In the case of bacterial sequencing one would expect >95% successful alignment, but when sequencing a mamallian genome (with many repetitive areas) it may be normal to have as low as 70-80% alignment success. RNA-Seq sequences regions that are usually well preserved, and thus alignment rates should be usually high. 

There can be several reasons why the alignment rate is low: the reads may not have been properly quality filtered or may contain artefactual sequence (such as adaptors and polyA tails); there may be contaminations; an inappropriate reference genome may have been used for alignment. The first of these reasons should have been detected before alingment, by looking at the raw data using FastQC and using the appropriate tools. It can be hard to find out if there were contaminations, unless we have an idea of the possible contaminants. An obvious one is Human, which we can check if we obtain a significant number of alignments against it. Finally, if we didn't use the proper genome, but there is no closer genome available, then there is not much that can be done, except perhaps trying to change parameters in the alignment software to allow for more mismatches (although this may cause biases and an increase in wrong alignments).

Another measure that can be used is the percentage of reads with duplicates (aligning exactly to the same place in the genome). Usually, duplication levels higher than 20% are not a good sign (they're a sign of low input DNA and PCR artifacts) but again, depends on what you are sequencing and how much. In RNA-Seq it is common to have a small set of genes highly expressed, leading to the presence of duplicates. The histogram of number of duplicates per read will often look bimodal, with most reads being unique and a small subset highly present (mostly from highly expressed genes). Unfortunately it is hard to distinguish PCR artifacts from highly expressed genes. When looking in IGV, PCR artifacts can be easily detected by an uneven coverage of the gene. To be safer, one can remove duplicates, but this is not usually done, since a lot of valid information may be lost.

QUESTION: Why duplication rates are frequently high in RNA-Seq? 

**Task**: In Galaxy, check the percentage of aligned reads in the alignments generated previously with sample datasets. Compare paired-end versus single-end, before and after trimming. Also compare the effect of changing the genome of reference.

Finally, there are reports specific for RNA-Seq which depend on gene annotation. One report indicates how well the genes are covered by sequence, which provides a good indication of RNA integrity. Finally, one can also check how well the alignments match the known annotation. The presence of a lot of alignments outside annotated genes can mean several things: annotation is not correct (eg. if you're working with a non-model organism); there can be DNA contamination; presence of immature RNA. Qualimap and [RSeqC](http://rseqc.sourceforge.net/) provide a set of tools to produce RNA-Seq specific reports. 

**Task**: Produce Qualimap (outside Galaxy) and RSseQC (in Galaxy) reports for the alignments you generated with your complete datasets. Some RSeqQC reports may take some time, so take care to run only one at a time during the day in Galaxy (similar to the alignments).

# Learning Outcome 7: Generate tables of counts using the alignment and a reference gene annotation

## LO 7.1 - The process of generating gene counts from genome aligments

To perform differential expression analysis we need to count, for each sample, how many times a different transcript/gene is read. If we align directly against the transcriptome, we just need to count the number of alignments per gene/transcript. However, if there are many alternative transcripts, aligning will become difficult. One solution may be to use just one representative transcript, or the union of all transcripts to represent the gene, although this also has issues.

What is most often done is to align against the genome, and compare the alignments (SAM/BAM) against the gene annotation (as GTF or BED). We could consider that a read counts to a gene if it overlaps with any part of the gene, but in large mammalian genomes, genes can have large introns, and it is not rare that genes overlap with each other. Moreover, the presence of DNA contamination and immature RNAs may also influence the counts.

Thus, it is usually preferable that a read will count for a gene only if it overlaps to at least some part corresponding to a valid mRNA transcribed from that gene. Then, if we have strand information, we should use it to resolve other possible ambiguities. But there are stil other factors to take in consideration. What to do if a read maps equally well to multiple genome regions? This will now depends a bit on the behavior on the alignment software. Usually, these cases are marked as having a low mapping quality, so we can simply ignore them by excluding alignments with a low mapping quality. But byt ignoring these cases we're losing information, and in the case of large genomes with a lot of large duplicated regions, this can be problematic. Again, if we want to use this information, we need to take into consideration what the aligner software will do. For example, bwa randomly attributes a read to one of the sites, while hisat outputs all alignmens (up to a given limit of k equally good ones). Some counting tools will actually use the information that a read aligns to different places to estimate the likelihood that a read belongs to one or the other, depending on the local (unique) coverage. This is in fact the type of approach Salmon uses to attribute reads to transcripts. Salmon does not output an exact number of reads per transcript, but the sum of the likelihoods of reads belonging to it (eg. a read may have 60% likelihood of belonging to a transcript, and thus will count not as 1, but as 0.6).

Finally, how to avoid pcr artifacts? To be as safe as possible, we would remove duplicates to avoid pcr artifacts, and this frequently needs to be done before the counting process. Nonetheless, given that duplicates can be frequent in RNA-Seq, usually we do not remove them. Assuming that pcr artifacts occurr randomly, then we should not have the same artifact in different biological replicates. In any case, for genes that are very important to use, we should always also visually check the alignments using software such as IGV.


## LO 7.2 - Use tools such as htseq-counts to generate table of gene counts

A popular tool to generate these counts from SAM/BAM alignments and GFF/GTF gene annotations is [htseq-count](http://www-huber.embl.de/HTSeq). Its default behavior is to generate counts at the gene level. It assigns a read to a gene if it unambiguously overlaps at least one part of a cDNA produced by the gene. It ignores reads mapping equally well to multiple positions by requiring by default a minimum mapping quality. By default it assumes stranded libraries, so we need to explicitly mention unstranded. 

**Task**: Use htseq-counts and Qualimap counts with the guilgur data (also use the sample gene annotations). Use the provided alignments (you may need to transform BAM to SAM to use with htseq-counts), and also try with your own alignments (obtained from the raw reads against a recent Drosophila genome). Try different parameters (namely related to strand and multiple mappings) and see the differences. Also check the effect of using an incorrect GTF file (not matching the correct genome version). 

**Task**: Run a Salmon "alignment" of the guilgur data against the sample transcriptome. Notice that no SAM/BAM is generated. Also compare results with the ones generated by htseq-counts. 

**Task**: Run htseq-count (and/or Qualimap counts) and salmon for your complete dataset. Similarly to the alignments, only run one at a time during the course and if necessariy leave them all running overnight.


# Learning Outcome 8: Generate lists of differentially expressed genes, at least for a simple pairwise comparison

## LO 8.1 - Using the R packages edgeR and DESeq2 to produce a pairwise differential expression analysis

The analysis methods currently most commonly used to perform RNA-Seq differential gene expression analysis start from non-normalized raw read counts like what we obtained previously (note that Salmon is not exactly like this). Given that sequencing data is based on discrete counts, most of these popular methods are based on derivations of the binomial distribution. Similarly to microarrays, there are many available tools to perform these analysis using the R language (such as [edger](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)).

The first thing that is done by edgeR or DESeq2 is to normalize the table of counts. Given that you can have different numbers of reads for each sample, you need to account for these differences. Nonetheless, this is not enough. You need to be aware that sequencing is a sampling experiment. If in one sample there is a gene or set of genes very highly expressed, then the amount of reads available for the rest of the genes is lower than for other samples where the highly expressed gene is not so highly expressed. So, many genes that may not be differentially expressed will appear so just because of this sampling artifact. 

If we assume that most genes are not differentially expressed, then the ratio of counts between samples should be the same for most genes. Under this assumption, we can obtain the mean counts for all genes, calculate the ratio against this reference mean for each gene in each sample, and then take the median of all ratios within one sample (to avoid the outliers). This median is then the sample specific normalization factor. This is the process that DESeq applies. EdgeR applies a similar, although more sophisticated, approach (trimmed mean of M-values, or TMM in short). 

To calculate differentially expressed genes, we need to take into consideration how much (normalized) counts vary between the different samples. This variation is clearly gene dependent, since highly expressed genes vary more in terms of absolute value, and low expressed genes vary more in terms of % of gene expression (fold change). If one only looks at fold change without taking variation into account, we’re more likely to have low expressed genes as differentially expressed. Therefore, we need to accurately estimate variation per gene, but we usually do not have enough replicates to do this on a gene by gene basis. One alternative that is used by edgeR and DESeq2 is to bin genes with similar expression and fit a curve.

We then test each gene for differential expression, and we obtain a probability for the test. Since we test thousands of genes, some genes may get good p-values just by chance. One way of avoiding this is by multiplying the p-value by the number of tests (a method called Bonferroni correction). This is nonetheless too strict and we usually end up not having anything differentially expressed. Other methods . We will look into more detail on this when we discuss functional enrichment analysis. In the end of a DESeq2 or edgeR analysis, instead of looking at the p-value, we should rather look at the corrected p-value (FDR, or qvalue) for significance. Finally, another way of minimizing the number of tests is to filter out the genes that have very low expression in all samples. 


## LO 8.2 - Interpretation and visualization of results

**Task**: In Galaxy, use DESeq2 with the htseq-count results you obtained previously for the guilgur data. Perform a simple parwise Wild-Type versus Mutant comparison with two replicates each. Look at the differentially expressed genes.

Even before interpreting the results of the differential expression analysis, we should have an idea of how the samples compare to each other. For this, we can look at plots such as Principal Coordinate Analysis (PCoA) or Multi-Dimensional Scaling (MDS). The way each software implements these plots vary a bit, but in short, they provide some evidence for the relatedness of the samples. Another common plot shows the hierarchical clustering of samples by displaying a heatmap with a matrix of distances between the samples. By looking at these plots we can detect outliers that may need to be removed from the analysis, or possible batch effects that we need to control.

DESeq2 and edgeR also plot the estimates of the biological coefficient of variation (BCV) plot, which depicts the sample variation of genes accordinf to their expression, and illustrates the variation correction the software performed, as we discussed in the previous section. Finally, another type of common plot is the MA or vulcano plot, which displays the average normalized expression of genes and their log fold change between the groups being compared. On top of these graphs it is common to signal the genes that were detected as differentially expressed. 

**Task**: In Galaxy, use DESeq2 to perform a pairwise comparison with the htseq-count results you obtained for your complete dataset. Look at the set of plots produced by the software.

Unfortunately, Galaxy does not produce gene-centered plots, and for those we may need to go to other software such as R. Nonetheless, the Galaxy tools output tables with normalized values that can be used for plotting in any type of software.


## LO 8.3 - Use more complex settings: Generalized Linear Models

So far, we just considered the simple case of pairwise comparison, where all samples are independent. But we may have cases where the samples are not independent. For example, in case of cancer, it is common (and desirable) to have tumor tissue and normal tissue for the same individual. In this case, we have paired information that needs to be taken into account in the test. There can also be other variables (eg. samples were prepared in different batches) that may confound the differential expression analysis. 

The [edgeR manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) contains several examples that explore these issues. In the pairwise case, the statistical methods are comparable to a t-test or a fisher exact test. Generalized Linear Models (GLM) allow to include information from several variables simultaneously. The simple pairwise case can also be considered as a  GLM, although the statistical methods applied for the test are different than in the "classic" pairwise model. 

In a first example, we have a classic paired test, with tumor and normal samples for one same case. In a second example, we have a treatment, but the samples were obtained in three different moments in time (batches), and this may influence the result.

**Task**: In Galaxy, use edgeR with the files provided for example1 and example2. In both cases, include paired/block information. Compare the difference with and without this information.

**Task**: In Galaxy, use edgeR to perform a pairwise comparison with the htseq-count results you obtained for your complete dataset. You will need to transform the htseq-count results into something edgeR can use.
		
The tools available in Galaxy are limited in terms of the ability to express more complex experimental designs. For this, we need to go to R and explore all the flexibility that it allows.
		
**Task**: Try the example1 and example2 using edgeR in R. For this, use the R script provided.

The final example we will explore contains several factors, and one of the factors have 3 different possible values. This introduces many possibilities of experimental questions to test. We just need to decide which ones are relevant biological questions.

**Task**: Try the example4 using edgeR in R. For this, use the R script provided.

**Task**: Based on what you did before, prepare a table of non-normalized counts for your complete dataset and analyse it using edgeR in R. 

# Learning Outcome 9 - Perform simple functional enrichment analysis and understand the concepts involved

## LO 9.1 - How to extract meaning from a list of genes

A list of genes of “interest” produced by an ‘omics experiment (e.g., RNAseq, microarrays, proteomics, etc) is essentially meaningless: gene identifiers are opaque, and we’re generally interested in understanding phenomena at the cellular and/or organismal level, rather than the gene level. To do this, we must abstract from the genes to their functions, or whatever other aspect we’re interested in studying (e.g., the chromosome location of the genes, the transcription regulation networks, etc).
In order to abstract to the functional level, we need functional descriptions of the genes. Furthermore, we need these descriptions to be consistent, i.e., we need all functional aspects to be described in the same manner for all genes that have those aspects − otherwise, we would be unable to integrate our gene set at the functional level. In short, we need genes to be annotated using a functional classification scheme.

There are several such schemes available, which cover different aspects and/or levels of protein function. For instance, the Enzyme Commission (EC) classification covers individual enzymatic function, whereas KEGG covers metabolic pathways, which is a different aspect (or view) of the same phenomena. However, there is only one classification scheme that covers a spectrum of gene function that is both wide and deep enough to analyze an ‘omics set as a whole: the Gene Ontology (GO).

GO is divided into three major functional aspects: molecular function, which covers individual gene functions; biological process, which covers how gene functions integrate into cellular and/or organismal processes; and cellular component, which covers where gene functions take place. Each of these aspects is organized as a directed acyclic graph, which is essentially a relaxed hierarchy with multi-parenting. In addition to subclass (‘is a’) relation, GO includes other relations such as ‘part of’, ‘occurs in’, and ‘regulates’. While the three aspects of GO are ‘is a’ orthogonal, molecular functions can be ‘part of’ biological processes, and both can ‘occur’ in cellular components.

GO is also available in the form of GO slims, which are ‘trimmed’ versions of the ontology where the specific fine grained terms have been removed and only broader terms are present. These usually cover the whole breadth of GO, albeit slims for particular species may exclude sections that are not applicable to that species. GO slims are useful for giving an overview of the GO annotations of a genome or a large collection of genes, when a broad classification is sufficient. However, they offer no advantage other than simplicity − whatever conclusion you derive using a GO slim would also be derived by using the whole ontology. They should not be used when a deeper classification is desired.

Genes can be (directly) annotated to multiple GO terms, even within the same aspect. Furthermore, according to the true path rule, a gene annotated to a GO term is implicitly annotated to all ancestors of that term. For instance, a gene annotated with ‘cytochrome c oxidase activity’ is inherently annotated with ‘catalytic activity’, ‘electron carrier activity’, ‘transporter activity’, and all other GO terms in the path between them.

GO annotations of genes are available, on an individual basis, in most genome databases, as well as in dedicate GO browsers such as AmiGO (http://amigo.geneontology.org) and QuickGO (https://www.ebi.ac.uk/QuickGO). They can also be downloaded on a genome-wide scale from GO’s annotation repository (http://www.geneontology.org/page/download-annotations) or BioMart (http://www.ensembl.org/biomart).
Viewing the annotations of your gene set on an individual gene basis is unfeasible and insufficient: there are too many genes to analyze manually and integrate, and even if you could, this doesn’t tell you how significant the patterns you find are.

**Task**: Go to BioMart through Galaxy (Galaxy > Get Data > BioMart) and get the GO annotations for the mouse genome in tsv (Gene Stable ID; GO term accession). Download the latest version of GO from http://geneontology.org/ontology/go.obo and upload it into Galaxy.       


## LO 9.2 - Understand the concept of functional enrichment analysis, and the statistics involved

Functional enrichment analysis is the application of Fisher’s exact test to measure the statistical significance of the observed frequency of each functional annotation in a gene set. The test relies on computing the probability of said frequencies arising by chance, using the hypergeometric distribution.

The statistical problem is essentially the same as determining the probability of getting x black balls in a sample of n balls, drawn without replacement from a bag with X black balls and a total of N balls. In our case the balls are genes, and a color is a functional annotation.

To compute the probability, we need to know the following parameters:
- The sample frequency: number of genes in the set annotated with the term
- The sample size: total number of genes (with any annotation) in the set
- The population frequency: number of genes in the population annotated with the term
- The population size: total number of genes (with any annotation) in the population

The population should be the total set of genes involved in the study: all expressed genes in the case of RNAseq, genes contained in the microarray in the case of microarray studies, etc. You should not use the whole genome as the population unless your dataset actually spans the whole genome. The study set can be defined in different ways: all genes with statistically significant differences in expression, only those with increased expression, or only those with decreased expression. It may make sense to perform enrichment analysis with all three study set options, as each gives you a different insight into your dataset. Note that you can also consider different thresholds for statistical significance, and different thresholds for expression increase/decrease, which may affect the enrichment analysis results.

When you use a GO enrichment analysis tool, you don’t need to define the sample or population frequencies. Instead, either you provide or the tool accesses the GO annotation set for your organism,  which is used to compute the frequencies for both sets. The tool should exclude from both sets the genes that don’t have GO annotations (of the aspect you are considering) − we don’t know their function at all, so we shouldn’t consider them as either positive or negative for any given functional aspect. However, not all tools do this, and it affects the results they produce.

Biological process is typically the most interesting aspect for enrichment analysis, but it may also be interesting to analyze the molecular function or cellular component aspect in particular studies − namely for validation. For instance, in a proteomics study where you sampled membrane proteins, you should do enrichment analysis with cellular component, to verify that “cellular membrane” is enriched.

Most of our statistical tests − including Fisher’s exact test − rely on controlling type I errors. When we accept an event/observation as significant because it has a p-value of 0.001, we are accepting that statistically, one time in a thousand, we’ll be wrong − the event/observation in question will be the product of chance alone. This is a problem when we perform multiple related tests, as the chance of getting a statistical “extreme” in at least one of them will be greater the more tests we perform. Because GO enrichment analysis relies on performing hundreds (or sometimes thousands) of Fisher’s tests, we must correct the statistics for multiple testing.

There are two families of multiple test corrections: the family-wise error rate (FWER) and the false discovery rate (FDR). In the former, we control the probability of making at least one false discovery, which is a conservative but “safe” approach that produces corrected p-values. In the latter, we control the ratio of false discoveries, which is a more powerful but less “safe” approach. It produces q-values, which indicate the ratio of false discoveries you are accepting if you reject the null hypothesis.

QUESTION: Why do we need multiple test corrections? What is the difference between a p-value, a corrected p-value, and a q-value?

Because GO is deeply hierarchic, performing enrichment analysis across GO requires a high number of tests, but not that many of them are independent. Thus, multiple test correction methods overestimate the error likelihood/rate. One way to reduce this effect is to not make redundant tests. For instance, if the frequencies of “protein binding” and its parent “binding” in your study set are the same, testing “binding” is redundant − the test can only be positive if the test of “protein binding” is positive, and the latter is more informative than the former.

There are several GO enrichment analysis tools available, for instance:
- Webtools: GOrilla (http://cbl-gorilla.cs.technion.ac.il/), GO’s own tool (http://www.geneontology.org/page/go-enrichment-analysis)
- Galaxy/stand-alone tools: GOEnrichment [IGC]; Ontologizer [test toolbox]
- R tools: gsea, GOstats, topGO

**Task**: Run an enrichment analysis test on GOrilla. Use the FEA_dataset1 from the Git repository, which contains the overexpressed genes from the Drosophila melanogaster dataset with 300 random genes differentially expressed. Choose “Two unranked lists of genes” as the running mode; paste or upload the study set into “target set” and the population set into “background set”. Use the biological process ontology, then repeat the analysis for the molecular function ontology.
Are there significantly enriched terms at 0.001 significance without multiple test corrections? And with the correction?


## LO 9.3 - Interpreting the results of functional enrichment analysis

It is essential to keep in mind that **statistically significant does not mean biologically meaningful**.

On the one hand, we can have functional enrichment of functional aspects that are too broad to derive any relevant conclusion, or that appear to be unrelated to the experiment in question. You should look at these with a critical eye − there may be some underlying meaning that is not readily apparent, even in the case of very generic terms such as “protein binding” or “single organism process”. In general, though, we’re more interested in functional aspects that are more specific.

On the other hand, aspects that are too specific may not be very interesting. In the extreme case of a biological process associated with a single gene in a given organism, if that gene appears in the study set, it is likely to be statistically enriched (if the study set is relatively small in comparison with the population set), but that doesn’t give us any insight into the study set as a whole.
In general, we’re interested in GO terms that are sufficiently generic to integrate a significant part of our dataset, but sufficiently specific to give us some conclusive insights.

Because of GO’s hierarchical structure, we may get related enriched terms with different levels of specificity, and we should consider them together as a cluster when drawing conclusions. These clusters may not be readily apparent from a results table, but are easy to detect in a graph view of the results (albeit graph views are not always easy to analyze due to the large size of GO).

It is also essential to consider that sporadic outliers may occur, despite multiple test corrections. Keep in mind that we’re making a statistical test (of enrichment) on top of another (of differential expression) which in turn is preceded by a statistical normalization. Even though we’re comfortable with the assumptions and p-values in each individual step, the likelihood of error propagates across the steps.

You should also keep in mind that enrichment analysis is qualitative, rather than quantitative: you are treating genes as either “on” or “off” (be “on” differentially expressed, overexpressed, or underexpressed) and consequently only assessing which functional aspects are statistically affected, rather than by how much they are affected.

**Task**: Run an enrichment analysis test on the GOEnrichment tool in Galaxy, using the FEA_dataset2 in the git repository. Use the go.obo and Mouse annotation file you got from BioMart earlier, as well as the dataset files. Run the program with the default options, then analyze the results tables.

**Task**: Run another enrichment analysis test on the GOEnrichment tool in Galaxy, using the FEA_dataset3 in the git repository. This time you will have to process the dataset to generate the lists of overexpressed and underexpressed genes. Hint: paste the dataset into a spreadsheet, and sort it and manipulate it there, then copy your up and down study sets, as well as your population set into text files. Analyze the results tables, then download the graph files and open them in yED.
