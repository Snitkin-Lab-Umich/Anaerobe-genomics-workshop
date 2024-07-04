Anaerobe-genomics-workshop
============================

Overview of workshop
--------------------
Over the past decade microbial genomics has made advances as the field has coallesced around best practices for standard analyses. Much of the standard workflows are comprised of open source software that is best (or exclusively) run on a command line (i.e. Linux) environment, ideally with a high-performance compute system (i.e. cluster). In this short session it is not feasible to learn how to navigate this environment and run all the tools neccesary for a genomic anlaysis (I teach a full semester course dedicated to this). So, instead we will walk through the steps in performing a genomic anlaysis, discuss how to interpret outputs and bring some data into R to get some experience working with processed data.

Overview of analysis
--------------------
We will be working through a data set that I published a few years back, which is associated with the following [paper](https://pubmed.ncbi.nlm.nih.gov/29167391/). The objective of this paper was to use whole-genome sequencig to understand the pathways of transmission in the context of a regional outbreak of carbapenem-resistant _Klebsiella pneumoniae_. This outbreak included 40 patients who were exposed to 26 different healthcare facilities. An epidemiologic investigation of this outbreak supported regional spread occuring via the movement of colonized and infected patients between different facilities. Our objective in this analysis was to provide further support for the hypothesis that patient transfer drives regional spread of antibiotic resistant organisms, and to see if genomic analysis could provide additional insight into the roles played by individual facilities in mediating regional spread, as well as the pathways by which a strain moved from one facility to another. 

Data generation
---------------
We performed whole-genome sequencing on isolates from outbreak patients. Steps included:

1. Plating patient samples on selective media with carbapenem antibiotics
2. Picking a single colony from culture plates to proceed with a single "clone" (i.e. single genotype) from the bacterial population in the patient
3. Innocuating the colony into broth culture to create additional material
4. Performing DNA extraction
5. Generating Illumina sequencing libraries that have adaptor sequences and barcodes
6. Sequencing isolates on Illumina MiSeq instruments

At many Universities there are dedicated sequencing cores that can perform steps 5 and 6 for you. In addition, there are commercial options available, whereby you can send bacterial pellets and have them not only sequenced, but also have data QC and basic analyses performed. We have done this using [plasmidsaurus](https://www.plasmidsaurus.com/), and been pleased with the results.

Data quality control
--------------------
A critical first step upon receiving your data is to make sure it is of sufficient quality, and take some additional steps to make sure it is ready for analysis. Things to take into account during data QC include:

1. Did you generate sufficient data for your analysis?
2. Is the sequence data of high quality (i.e. individual base calls made with high confidence)?
3. Is your sequence pure (i.e. no contaminants) and the expected organism?

There are different tools available to aid with these assessments that are standardly used in the field. 

### Overall sequence quality assessment
For evaluating the amount and quality of data, it is common to use a tool called [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC takes as input the fastq files that comprise the millions of short (100-200bp) reads that are generated for each genome, and provides different assessments of sequence quality as well as provide evidence for contamination.

### Read trimming
Oftentimes your data is overall good quality, but there may be a few low quality sequences, or residual adpaptor sequence that needs to be trimmed off. For this reason, even if FastQC indicates your data is of good quality, you should run a tool called [trimmomatic](https://github.com/timflutre/trimmomatic). Trimmomatic will: i) trim off adaptor sequences, ii) remove low quality reads and iii) trim off low quality parts of otherwise high quality reads.

### Genome assembly and assembly evaluation
After you have assessed the overall quality of your data and trimmed it, there are some additional tools to make sure you sequenced the right thing and that there is not any contamination that would cause problems in downstream analysis. To proceed with these assessments, a first step is to create a genome assembly. A genome assembler will take your short reads, and stitch them together to form larger contiguous sequences (i.e. contigs), which are suitable for downstream analyses of genome structure and function. The most common genome assembler for bacterial genomes is call Spades: here is a nice tool that makes it more user freindly called [Shovil](https://github.com/tseemann/shovill). After you construct an assembly, you should evaluate it's quality. Three things to pay attention to are the number of contigs, the overall size of the assembly and the typical size of contigs. A large number of contigs (i.e. > 500) is an indication of potential contamination. For the genome size, you should verify that it matches what is typical for the species you are intending to sequence. The typical size of contigs (i.e. N50) should be at least 20 Kb, as smaller than that is a sign of contamination. These assembly metrics can be generated with a tool called [quast](https://github.com/ablab/quast).

### Species determination
Once you have a genome assembly, you can use a tool called [skanni](https://github.com/bluenote-1577/skani) to compare to a database of genomes to assess the most likely species.

### Genome completeness
Lastly, you can use a tool called [CheckM](https://github.com/Ecogenomics/CheckM) to evaluate genome completeness. This tool works by verifying that your assembly has one and exactly one of known single copy genes. Not having all of these is a sign of insufficient sequencing, and having multiple copies is an indication of contamination.

Typing genomes using multi-locus sequence typing (MLST)
-------------------------------------------------------


Annotating antimicrobial resistance (AMR) genes and mutations with AMRfinderPlus
--------------------------------------------------------------------------------


Downloading public data to put local genomes in context
-------------------------------------------------------


Identifying variants in outbreak genomes and constructing a phylogeny
---------------------------------------------------------------------


Overlaying meta-data on your whole-genome phylogeny to understand outbreak
--------------------------------------------------------------------------


