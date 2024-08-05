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

## Copy required files to your directory
### copy fastq files
```
cd data
cp -R /mnt/proj/omicss24/alignment/data/* .
```
### copy reference genome
```
cd ../ref_genome
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
../data/SRR11881059_1.fastq.gz \
../data/SRR11881059_2.fastq.gz > aln.sam
```

## Sort the sam file and convert it to bam
```
samtools sort aln.sam -@2 -o aln_sorted.bam
```

## Create an index of a bam file
```
samtools index aln_sorted.bam
```
