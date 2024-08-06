# Alignment

## Connect to the server
### Please replace username with your uid provided in Slack
```
ssh username@comp1.genomic.abi.am
```
### Go to the working directory
```cd /mnt/user/username```

## Create a folder for the practice and enter the folder
```
mkdir alignment_task
cd alignment_task
```

## Create project folders
```
mkdir ref_genome # directory to store reference genome and its index files
mkdir data # directory to store fastq files
mkdir aln_res # directory to store alignment results
```

## Download and process required files for alignment
### Download fastq files with SRA toolkit
```
cd data
fastq-dump --gzip --skip-technical --split-files --clip SRR11881059
```
### Perform fastqc on downloaded files
```
fastqc *
```
### Trim illumina adapters
```
mkdir trimmed
cd trimmed
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o SRR11881059_1_trimmed.fastq.gz -p SRR11881059_2_trimmed.fastq.gz ../SRR11881059_1.fastq.gz ../SRR11881059_2.fastq.gz
```
### Run fastqc again
```
fastqc *
```

### copy reference genome
```
cd ../../ref_genome
cp /mnt/proj/omicss24/alignment/ref_genome/GCA_000146045.2_R64_genomic.fna .
```

## index reference genome
```
cd ../ref_genome
bwa index GCA_000146045.2_R64_genomic.fna
cd ..
```

## Perform alignment with bwa mem
```
cd aln_res
bwa mem ../ref_genome/GCA_000146045.2_R64_genomic.fna -t 2 \
../data/trimmed/SRR11881059_1_trimmed.fastq.gz \
../data/trimmed/SRR11881059_2_trimmed.fastq.gz > aln.sam
```

## Sort the sam file and convert it to bam
```
samtools sort aln.sam -@2 -o aln_sorted.bam
```

## Create an index of a bam file
```
samtools index aln_sorted.bam
```
