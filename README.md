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

#### Overall sequence quality assessment
For evaluating the amount and quality of data, it is common to use a tool called [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC takes as input the fastq files that comprise the millions of short (100-200bp) reads that are generated for each genome, and provides different assessments of sequence quality as well as provide evidence for contamination.

#### Read trimming
Oftentimes your data is overall good quality, but there may be a few low quality sequences, or residual adpaptor sequence that needs to be trimmed off. For this reason, even if FastQC indicates your data is of good quality, you should run a tool called [trimmomatic](https://github.com/timflutre/trimmomatic). Trimmomatic will: i) trim off adaptor sequences, ii) remove low quality reads and iii) trim off low quality parts of otherwise high quality reads.

#### Genome assembly and assembly evaluation
After you have assessed the overall quality of your data and trimmed it, there are some additional tools to make sure you sequenced the right thing and that there is not any contamination that would cause problems in downstream analysis. To proceed with these assessments, a first step is to create a genome assembly. A genome assembler will take your short reads, and stitch them together to form larger contiguous sequences (i.e. contigs), which are suitable for downstream analyses of genome structure and function. The most common genome assembler for bacterial genomes is call Spades: here is a nice tool that makes it more user freindly called [Shovil](https://github.com/tseemann/shovill). After you construct an assembly, you should evaluate it's quality. Three things to pay attention to are the number of contigs, the overall size of the assembly and the typical size of contigs. A large number of contigs (i.e. > 500) is an indication of potential contamination. For the genome size, you should verify that it matches what is typical for the species you are intending to sequence. The typical size of contigs (i.e. N50) should be at least 20 Kb, as smaller than that is a sign of contamination. These assembly metrics can be generated with a tool called [quast](https://github.com/ablab/quast).

#### Species determination
Once you have a genome assembly, you can use a tool called [skanni](https://github.com/bluenote-1577/skani) to compare to a database of genomes to assess the most likely species.

#### Genome completeness
Lastly, you can use a tool called [CheckM](https://github.com/Ecogenomics/CheckM) to evaluate genome completeness. This tool works by verifying that your assembly has one and exactly one of known single copy genes. Not having all of these is a sign of insufficient sequencing, and having multiple copies is an indication of contamination.

Typing genomes using multi-locus sequence typing (MLST)
-------------------------------------------------------
Before performing an in depth genomic analysis to characterize the relatedness among your isolate genomes, it is often helpful to perform a higher level assessment. In the context of molecular epidemiology, a common approach for characterizing and cataloging genetic diversity of a strain collection is multi-locus sequence typing (MLST). MLST works by determining the sequence of a predetermined set of 6-7 conserved genes, and assigning a type based on the sequence of these genes. In the context of an outbreak investigation, this can be an informative first step. If all of your isolates are part of a single outbreak, then you'd expect them to be all of the same ST. Similarly, if there are some isolates that are different from the dominant ST, then they are unlikely to be related by clonal transmission (i.e. not part of the outbreak). A second advantage of performing MLST as a first step, is it allows you to put your genomes into buckets for downstream high-resolution genomic anlayses, where you can identify genetic variation among groups of genomes that belong to the same ST, which is usually of primary interest.

When sequencing the whole genome, you are also getting the sequence of the subset of MLST genes, allowing you to infer ST from genomic data. A common tool for extracting MLST from genomes is called [mlst](https://github.com/tseemann/mlst).

Applying this tool to our CRKP outbreak genomes reveals that they are all ST 258, which is the most common strain of CRKP in the United States. This supports the regional isolates indeed being linked by transmission.

Annotating antimicrobial resistance (AMR) genes and mutations with AMRfinderPlus
--------------------------------------------------------------------------------
When investigating bacterial pathogens it is often of interest to catalog antibiotic resistance and virulence factors present in genomes. A common tool to perform this task is [AMRfinderplus](https://github.com/ncbi/amr). AMR finder plus compares genes from your genome to a curated database of resistance/virulence genes, as well as resistance mutations for several epidemiologically important organisms.

Applying AMRfinderplus to our outbreak genomes shows that they have many antibiotic resistance genes, as one might expect for a healthcare associated lineage like ST258. Of greatest interest is the presence of the _Klebsiella pneumoniae_ carbapenemase (KPC) enzyme, which confers carbapenem resistance and makes CRKP infections difficult to treat.

In addition to annotating AMR and virulence determinents, a standard part of any genomics pipeline is to perform broader genome annotation (i.e. detect all putative genes encoded in a genome). The most commonly used tool for this is [Prokka](https://github.com/tseemann/prokka), which will detect protein coding and RNA genes, and then compare predicted genes to databases of known genes to assign putative functions.

Downloading public data to put local genomes in context
-------------------------------------------------------
Oftentimes it useful to put genomes that you have sequenced in the context of publically available isolates. Example analyses of interest include:

1. Understanding where in the world geneticlly similar strains have been observed, to gain insight into the potential origin of your strians.
2. Determining whether there are unique characteritics of your strains as compared to what has been observed previously (i.e. novel resistance, virulence, etc.)

In the context of our outbreak investigation, we suspect that it is indeed a clonal outbreak based on all isolates being the same ST. However, we can garner further support by seeing how our isolates compare to public isolates. If our isolates indeed derive from clonal spread in the region, then we expect then to be more closely related to each other, than to isolates we see in public databases. There are several databases that house microbial genomes, including [PATRIC](https://www.bv-brc.org/) and [NCBI](https://www.ncbi.nlm.nih.gov/). Here we will get genomes from NCBI from a previously published paper characterizing the genetic diversity of [ST258](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0133727). When papers with genomic data are published, they are required to submit them to sequence databases. When submitting a genomic dataset to NCBI it is typically associated with a bioproject, which will be indicated in the original manuscript. We can then use this bioproject ID to pull the genomes of itnerest from NCBI.

To get a list of sample IDs associated with this bioproject we will use tools provided by NCBI. In particular, NCBI provides a suite of command line tools called Entrez Direct also known as E-utilities to access the metadata stored in its various databases. The three main tools of E-utilities are - esearch, esummary and xtract that lets you query any of the NCBI databases and extract the metadata associated with the query. The query can be anything(Bioproject, Biosample, SRA accession, Genbank assembly accession).

Here is the command that we used to extract metadata information for the above mentioned research study.

```
esearch -db sra -query PRJNA252957 | esummary | xtract -pattern DocumentSummary -element Experiment@acc,Run@acc,Platform@instrument_model,Sample@acc > PRJNA252957-info.tsv
```

We can then use the IDs present in this file to download the genomes of interest using the ncbi tool fasterq-dump. 

```
cut -f2 PRJNA252957-info.tsv | parallel fasterq-dump {}
```

Once we have downloaded a large dataset, we want a way to quickly describe the genetic relationships. Below we will perform a fine-grained analysis of the relationship among our outbreak genomes, where accuracy of genetic distances is critical. However, for this sort of contextual anlaysis, we often prioritize the size of datasets and the speed of analysis, versus the accuracy. To get approximate genetic relationships among our downloaded and outbreak genoems we can use a tool called [mashtree](https://github.com/lskatz/mashtree).

Mashtree creates a phylogenetic tree using the neighbor joining (NJ) algorithm. We can use the web tool [iTOL](https://itol.embl.de/) to visualize this tree, along with meta-data describing where each isolate came from (i.e. outbreak or public). You can see that all of our outbreak genomes group together on the tree, again providing strong support that all cases are linked by transmisison in the region.


Identifying variants in outbreak genomes and constructing a phylogeny
---------------------------------------------------------------------
One of the most common goals of in sequencing a microbial genome is to identify small genetic variants like nucleotide substitutions or small insertions/deletions (i.e. indels). Identifying these genetic variants in one or multiple genomes has several important downstream applications:

1. Phylogenetic analysis - The first step in a phylogenetic analysis is identification of single nucleotide variants (SNVs) across the set of genomes of interest. Essentially, the input to any phylogenetic tree building software is a variant alignment that indicates what nucleotide each genome has at a position that is variable in at least one of the input genomes. The resulting phylogenetic tree then groups genomes together based on shared evolution, which is inferred from shared variants.

2. Transmission analysis - One of the most common approaches to assess the confidence in a putative transmission linkage between two individuals is to count the number of variants between two individuals pathogen genomes. Thus, having accurate variants is essential to making correct transmission inferences.

3. Functional analysis - In previous sessions we identified different types of genetic variation including differences in gene content and antibiotic resistance variation. Small changes in the genome like those identified through variant calling pipelines can also have significant functional impacts on genome function.

A typical variant calling analysis involves:

**Read Mapping**: Mapping sequenced reads to the reference genome using a read mapper

**Variant calling**: Calling variants(differences) between the reference genome and our sample.

**Variant filtering**: Filtering out variant calls that are deemed low confidence based on user defined criteria.

**Variant annotation**: Annnotating these variants to learn about their their effect on proteins and other biological processes.

There are many pipelines that stitch all of these steps together, with among the most popular due to its accuracy and ease of use being [snippy](https://github.com/tseemann/snippy). Using the output of snippy you can:

1. Perform a phylogenetic analysis by feeding directly into a tool like [gubbins](https://github.com/nickjcroucher/gubbins). Gubbins performs recombination masking (i.e. identifies and masks putative horizontally acquired variants) and then performs maximum likelihood phylogenetic analysis.

2. View read alignments and variant calls using a graphical tool like [IGV](https://igv.org/).
   
3. Examine annotations of variants using functional annotation provided by snippy 

Overlaying meta-data on your whole-genome phylogeny to understand outbreak
--------------------------------------------------------------------------
Having generated a phylogeny from genomic variants, we next will visualize it with labels for the facility from which isolates came from. We could do this in iTOL, but instead I will show you how this can be done in R, which enables more customizable, reproducible and quicker analysis.


First let's plot a bar plot to see how many isolates we have from each facility.
```
#Load libraries
library(ape)
library(ggplot2)
library(ggtree)

#Read in meta-data
meta_data <- read.table('assignment_3/crkp_outbreak_metadata.txt',
                        sep = "\t",
                        header = T)

#Read in tree
tree <- read.tree('assignment_3/crkp_gubbins.tree')

#Read in variant alignment
aln <- read.dna('assignment_3/crkp_gubbins_var_aln.fasta',
                format = "fasta")


#Plot distribution of isolates across facilities
ggplot(data=meta_data, aes(x=Isolation_Source)) +
  geom_bar(stat="count") + 
  theme(axis.text.x = element_text(angle = 90))
```

Next let's use [ggtree](https://yulab-smu.top/treedata-book/) to plot the phylogeny with facility tip labels.

```
#Get isolates named by facility
facils <- structure(meta_data$Isolation_Source, names = meta_data$Run)

#Get vector of facilities for tips and nodes of tree
facils_tip <- c(facils[tree$tip.label], 
                rep(NA, Nnode(tree)))

#Plot tree with ggtree
ggtree(tree) + 
  geom_tippoint(aes(col = facils_tip)) + 
  scale_color_discrete() +
  labs(col = 'Facility')
```




