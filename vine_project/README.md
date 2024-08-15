# Grapevine Project

In this project, we will reanalyze data from [this publication](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1386225/full). The study focuses on identifying genetic markers in the grapevine genome associated with the resistance to downy mildew (DM) and powdery mildew (PM). DM and PM are fungal diseases that significantly affect grapevines, leading to substantial crop losses. Breeding grapevines for resistance to these diseases is crucial for sustainable viticulture. The researchers used Genome-Wide Association Analysis (GWAS) to identify Single Nucleotide Polymorphisms (SNPs) associated with resistance. The study successfully identified new loci linked to DM and PM resistance.

In this project, we will replicate the main steps of the study, which include:

- **Data Preprocessing** - Filtering SNPs and samples according to the provided pipeline.
- **PCA Analysis** - Conducting Principal Component Analysis (PCA) to explore genetic variation.
- **ADMIXTURE Analysis** - Performing ADMIXTURE analysis to determine population structure.
- **GWAS Analysis** - Conducting GWAS to identify SNPs associated with DM and PM resistance in grapevines.
- **Functional Annotation** - Annotating the SNPs significantly associated with resistance.
















- **PCA**

input: /mnt/proj/vine/shared_files/SSFinalProj/data/plink/genotype
output: pca plots

1.Perform  PCA analysis in the commad line using the tool plink, get eigenvec and eigenvalue files

2.Visualize the PCA results using R, colour the samples by subspecies



- **Performing quality check filtering of the PLINK files before doing GWAS**

Read more about the quality check steps here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/

input: /mnt/proj/vine/shared_files/SSFinalProj/data/plink/genotype
output: genotype_filtered PLINK files

Do the following using the plink tool (keep track of how many SNPs and samples remain after each filtering step).

1. Delete SNPs and samples with high levels of missingness:

Delete SNPs with missingness percentage over 2%.

Delete samples with missingness percentage over 2%.

2.Remove SNPs with a low MAF frequency. (A conventional MAF threshold for a regular GWAS is 0.05)

3. Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE):

By default the --hwe option in plink only filters for controls.
Therefore, we use two steps, first we use a stringent HWE threshold for controls (1e-6), followed by a less stringent threshold for the case data (1e-10).

4.  Remove sampless with a heterozygosity rate deviating more than 3 sd from the mean.

Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
Therefore, to generate a list of non-(highly)correlated SNPs,  prune the SNPs using the command --indep-pairwise.
The parameters  50 5 0.2  stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.

Afterwards, you can  generates a list of samples that deviate more than 3 standard deviations from the heterozygosity rate mean using R and remove them.



Filtering is done!


- **Convert the quality check passed plink file to VCF using the tool plink**

- **Convert the VCF file to hapmap format to use in GAPIT using the tool TASSEL**

- Delete SNPs on unknown chromsomes (chr0) and the depricated 606th row

- **Run GAPIT GWAS models (GLM, MLM,CMLM,SUPER,MLMM,FarmCPU,BLINK) using the GAPIT package in R**

- Read more about the GAPIT tool and models here: https://zzlab.net/GAPIT/gapit_help_document.pdf

- **Functioanl Annotation**
  1. Find the genes that are located 0.9 mb upstream and downstream the significantly associated SNPs from the annotation gff file.
  2. Do gene enrichement analysis using DAVID enrichment web tool: https://david.ncifcrf.gov/tools.jsp 
