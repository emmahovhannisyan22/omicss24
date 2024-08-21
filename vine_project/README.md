# Grapevine Project

In this project, we will reanalyze data from [this publication](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1386225/full). The study focuses on identifying genetic markers in the grapevine genome associated with the resistance to downy mildew (DM) and powdery mildew (PM). DM and PM are fungal diseases that significantly affect grapevines, leading to substantial crop losses. Breeding grapevines for resistance to these diseases is crucial for sustainable viticulture. The researchers used Genome-Wide Association Analysis (GWAS) to identify Single Nucleotide Polymorphisms (SNPs) associated with resistance. The study successfully identified new loci linked to DM and PM resistance.

In this project, we will replicate the main steps of the study, which include:

- **Data Preprocessing** - Filtering SNPs and samples according to the provided pipeline.
- **PCA Analysis** - Conducting Principal Component Analysis (PCA) to explore genetic variation.
- **ADMIXTURE Analysis** - Performing ADMIXTURE analysis to determine population structure.
- **GWAS Analysis** - Conducting GWAS to identify SNPs associated with DM and PM resistance in grapevines.
- **Functional Annotation** - Annotating the SNPs significantly associated with resistance.

## GWAS pipeline

In this project, we will use data from the publication. The genotype data is provided as a table in [Supplementary Table 2](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1386225/full#supplementary-material), and the phenotype information is stored in [Supplementary Table 3](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1386225/full#supplementary-material). Take a look at the raw data.
To make things easier, we have already transformed the provided matrix into the commonly used PLINK binary format, which will be used in the pipeline. The files are stored on the server in the following directory:

> /mnt/proj/omicss24/vine_project/data/plink

The plink files consists of 3 following seperate files.
- genotype.bed
- genotype.bim
- genotype.fam


Get familiar with the file format in the  [plink documentaion](https://www.cog-genomics.org/plink/1.9/formats). 

### Directory setup

During this project, you will work in your own directory. Let's set it up!

> cd /mnt/user/username \
> mkdir vine_project \
> cd vine_project/ \
> mkdir scr # directory for scripts \
> mkdir scr/log #directory for the log files 
> mkdir result # directory for the results \
> mkdir data/metadata \
> mkdir data # directory for data \
> cp -r /mnt/proj/omicss24/vine_project/data/plink/ data/ # copy input files \
> cp /mnt/proj/vine/shared_files/SSFinalProj/data/metadata/viticola_pheno.txt data/metadata \
> cp /mnt/proj/vine/shared_files/SSFinalProj/data/metadata/necator_pheno.txt data/metadata \
> cp /mnt/proj/vine/shared_files/SSFinalProj/data/metadata/subspecies_pheno.txt data/metadata

Now you are ready to go! 

### Performing quality check filtering of the PLINK files before doing GWAS

Read more about the quality check steps [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)

All the preprocessing steps should be done by **plink** tool. It is already installed on the server. Here is the [plink documentation](https://www.cog-genomics.org/plink/) 

input: **data/plink/genotype**

output final: **data/plink/genotype_filtered** PLINK files

Do the following using the plink tool (keep track of how many SNPs and samples remain after each filtering step). After each filtering step the **.log** file is generated, where you can track the progress. 

1. Delete SNPs and samples with high levels of missingness:

Delete SNPs with missingness percentage over 2%.

Delete samples with missingness percentage over 2%.

2. Remove SNPs with a low MAF frequency. (A conventional MAF threshold for a regular GWAS is 0.05)

3. Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE):

By default the --hwe option in plink only filters for controls.
Therefore, we use two steps, first we use a stringent HWE threshold for controls (1e-6), followed by a less stringent threshold for the case data (1e-10).

4.  Remove samples with a heterozygosity rate deviating more than 3 sd from the mean.

4.1 Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
Therefore, to generate a list of non-(highly)correlated SNPs,  prune the SNPs using the command --indep-pairwise.
The parameters  50 5 0.2  stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.
Note: the pruned dataset is only used for heterozygosity checks, not the further steps.

4.2 Computes observed and expected autosomal homozygous genotype counts for each sample (--het)

4.3 Using R, get heterozygosity rates of each sample. Note: 

$$
\text{Heterozygosity rate} = \frac{\text{O(HET)}}{\text{N(NM)}} = \frac{\text{N(NM)} - \text{O(HOM)}}{\text{N(NM)}}
$$

4.4 Using R, calculate mean and standard deviation of hetrozygosity rate. Find and exclude samples which deviate more than 3 sds from the mean.
​
 

Filtering is done!

## PCA 


input: data/plink/genotype
output: pca plots


1.Perform  PCA analysis in the commad line using the tool plink, get eigenvec and eigenvalue files

2.Visualize the PCA results using R, colour the samples by subspecies


## ADMIXTURE

1) Create vine_project/soft folder and
2) Download the ADMIXTURE software from here (https://dalexander.github.io/admixture/download.html)
3) Unzip the .tar.gz file, In this folder you can find also the Manual
4) Run ADMIXTURE for k= 1-10, you can write a for loop for this. The input file should be the filtered, pruned plink file (data/plink/genotype_maf005_hwe.bed)
5) Enable the cross-validation calculation with the -cv option
6) Save the log files in the scr/log directory by adding  " | tee ./scr/log/log${K}.out "
7) Explore the output files, you can read about them in the manual.
8) Grep the CV errors from the log files, save them in a separate file and plot them using R
9) For plotting the ADMIXTURE results you can use this code (https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r)



## Preparing data format for GWAS 

- **Convert the quality check passed plink file to VCF using the tool plink**

- **Convert the VCF file to hapmap format to use in GAPIT using the tool TASSEL**

- Delete SNPs on unknown chromsomes (chr0) and the depricated 606th row



## GWAS analysis

- **Run GAPIT GWAS models (GLM, MLM,CMLM,SUPER,MLMM,FarmCPU,BLINK) using the GAPIT package in R**

- Read more about the GAPIT tool and models here: https://zzlab.net/GAPIT/gapit_help_document.pdf

## Functional Annotation

  1. Find the genes that are located 0.9 MB( linkage disequilibrium range form the paper) upstream and downsteram the singificantly associated SNPs.
  2. Do gene enrichment analysis using [DAVID web tool](https://david.ncifcrf.gov/tools.jsp)

