# Overview

### Method

1. Align reads to genome and transcriptome using Star
2. Count genome mappings using HTSeq
3. Count transcript mappings using EdgeR and generate counts
4. Generate bigwig files for reads mapping to genome

# How to run the pipeline

1) Get the code

```
git clone https://github.com/elswob/rna-seq-pipe.git
```

2) Put all your FASTQ files in a fastq folder and name the folders with the sample names or how you want them.

3) Create a txt file called Conditions.txt with folder/sample names and conditions in a tab delimited file

```
Sample1	Condition1
Sample2 Condition1
Sample3 Condition2
Sample4 Condition2
Sample5 Condition3
Sample6 Condition3
```

4) Check setup - there should be three things in the directory:

1. A folder containing the fastq files in separate sample folders with names corresponding to those in conditions file
2. A Conditions.txt file as above
3. A symlink to the pipeline

5) To run the scripts just Run run.sh with the following arguements:  

1. Number of CPUs  
2. Species (Mouse or Human)  
3. Location of a tmp directory
4. Location of directory containing samples and fastq files  
5. Strandedness (reverse, forward or no)

```
./combined/Run.sh 5 Mouse /path/to/temp/folder/ fastq_folder/ reverse
```

