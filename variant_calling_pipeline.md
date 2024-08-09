# Variant Calling Pipeline with GATK

This pipeline leverages the [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/) to identify variants in genomic data. The pipeline covers preprocessing steps, variant discovery, and filtering to ensure accurate and reliable variant calls. An example using the human MEFV gene sequence for one sample is provided. This pipeline is based on the [GATK Best Practices Workflows](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).

# Directory Setup
In your home directory, copy the following folder with all data and scripts:

> cp -r /mnt/proj/omicss24/variant_calling ./

This folder contains the following directories and files:
- **data/fastq** - Folder containing FASTQ files of the Illumina paired-end sequences for one sample.
- **data/reference** - Folder containing the human reference genome in FASTA format, along with the BWA index.
- **scr/variant_calling_pipeline.sh** - Script containing all commands for variant calling.
- **soft** - All the necessary soft

# Alignment
Raw reads need to be aligned to the human genome. First, an index should be created for the reference genome. Since genome indexing is time-consuming, it has already been performed, so you do not need to do it again. The code for indexing is provided below for reference:

> bwa index data/reference/hg38.fa # Genome indexing

**Note:** The read group information is where you enter the metadata about your sample.

Create a directory for the alignment results:

> mkdir data/bam_sam  # Creating a directory for the alignment results

Align the raw reads using the BWA tool:

> bwa mem \
> -R '@RG\tID:sample1\tLB:sample1\tSM:sample1\tPL:ILLUMINA’ \
> data/reference/hg38.fa \
> data/fastq/HbFMF_1_MEFV.fastq  \
> data/fastq/HbFMF_2_MEFV.fastq  > data/bam_sam/HbFMF.sam

## Sorting SAM File and Marking Duplicates

Sorting the SAM file is a necessary step before marking duplicates. During the sequencing process, the same DNA molecules can be sequenced several times. The resulting duplicate reads are not informative and should not be counted as additional evidence for or against a putative variant. The duplicate-marking process (sometimes called "dedupping" in bioinformatics slang) identifies these reads so that the GATK tools know to ignore them.

As a result, you will get a sorted BAM file called **HbFMF.sort.markdup.bam**. This file contains the same content as the input file, except that any duplicate reads are marked. Here are the essential steps:
- **Sorting:** The input sequencing data is first sorted by genomic coordinates.
- **Identifying duplicates:** Duplicates are defined as reads that have the same start and end coordinates.
- **Marking duplicates:** Once duplicates are identified, one of the reads is marked as the primary or representative read, while the duplicates are marked as secondary. This marking process involves modifying the flags or tags associated with each read in the BAM file.
- **Optional duplicate removal**

Sort the SAM file:

> /mnt/proj/omicss24/variant_calling/soft/gatk-4.6.0.0/gatk SortSam \
> --INPUT data/bam_sam/HbFMF.sam \
> --OUTPUT data/bam_sam/HbFMF.sort.bam \
> --SORT_ORDER coordinate

Mark duplicates:

> /mnt/proj/omicss24/variant_calling/soft/gatk-4.6.0.0/gatk MarkDuplicates \
> -I data/bam_sam/HbFMF.sort.bam \
> -O data/bam_sam/HbFMF.sort.markdup.bam \
> -M data/bam_sam/HbFMF.sort.markdup_metrics.txt

Index the BAM file:

> samtools index data/bam_sam/HbFMF.sort.markdup.bam

Now you can download the **HbFMF.sort.markdup.bam** file and visualize it in IGV.

# Variant Calling

Create a directory for VCF files:

> mkdir data/vcf # Create directory for VCF files

Run variant calling:

> /mnt/proj/omicss24/variant_calling/soft/gatk-4.6.0.0/gatk -T HaplotypeCaller \
> -R data/reference/hg38.fa \
> -I data/bam_sam/HbFMF.sort.markdup.bam \
> -o data/vcf/raw_variants.vcf

# Variant Filtration

### Separating SNPs and INDELs
The result of variant calling is the **raw_variants.vcf** file, which includes both SNPs and small INDELs. However, in later steps of variant filtration, they should be separated. Create two separate VCF files for SNPs and INDELs:

For SNPs:

> /mnt/proj/omicss24/variant_calling/soft/gatk-4.6.0.0/gatk select-type SNP \
> -V data/vcf/raw_variants.vcf \
> -o data/vcf/snp_variants.vcf

For INDELs:

> /mnt/proj/omicss24/variant_calling/soft/gatk-4.6.0.0/gatk select-type INDEL \
> -V data/vcf/raw_variants.vcf \
> -o data/vcf/indel_variants.vcf

### Filtering Variants

In the default parameters, the most commonly used filters for SNPs and INDELs are as follows:

- **QD (Variant Confidence/Quality by Depth):** The Variant Confidence/Quality by Depth (QD) annotation in GATK is a metric used to estimate the confidence of a variant call. It quantifies the quality of a variant by normalizing it with respect to the depth of sequencing coverage at that position.
- **QUAL (Quality):** This field in the Variant Call Format (VCF) file represents the Phred-scaled quality score of a variant call.
- **SOR (Strand Odds Ratio):** SOR is an annotation that quantifies the strand bias observed in the reads supporting a variant call. It measures the imbalance in the distribution of reads across the forward and reverse strands at a variant site.
- **FS (FisherStrand):** FS is an annotation that quantifies strand bias using Fisher's exact test. It measures the probability of observing the distribution of reads supporting the variant call across the forward and reverse strands under the assumption of no strand bias.
- **MQ (Mapping Quality):** MQ is assigned to each variant based on the mapping qualities of the reads aligned to the variant position.
- **MQRankSum:** MQRankSum is calculated by ranking the mapping qualities of the reads that support each allele and then calculating the difference in ranks between the alternate and reference alleles.
- **ReadPosRankSum:** ReadPosRankSum is calculated by ranking the positions of the reads that support each allele and then calculating the difference in ranks between the alternate and reference alleles.

#### SNPs

Apply the following filters for SNPs:

> /mnt/proj/omicss24/variant_calling/soft/gatk-4.6.0.0/gatk -T VariantFiltration \
> -R data/reference/hg38.fa \
> -V data/vcf/snp_variants.vcf \
> -filterName "QD_filter" \
> -filter "QD<2.0" \
> -filterName "FS_filter" \
> -filter "FS>60.0" \
> -filterName "MQ_filter" \
> -filter "MQ<40.0" \
> -filterName "SOR_filter" \
> -filter "SOR>10.0" \
> -o filtered_snps.vcf

#### INDELs

Apply the following filters for INDELs:

> /mnt/proj/omicss24/variant_calling/soft/gatk-4.6.0.0/gatk -T VariantFiltration \
> -R data/reference/hg38.fa \
> -V data/vcf/indel_variants.vcf \
> -filterName "QD_filter" \
> -filter "QD<2.0" \
> -filterName "FS_filter" \
> -filter "FS>200.0" \
> -filterName "SOR_filter" \
> -filter "SOR>10.0" \
> -o filtered_indels.vcf

#### End of Pipeline ####
