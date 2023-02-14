# biol217 (Marine)Microbial Omics - from sample to function

Aim of the course:
- Practical and theoretical knowledge of a bioinformatic workflow
- Assembly and analyse of MAGs 
- Skills in the command line interface
- documentation of bioinformatic analyses

## Dataset
We worked with samples taken from a mesophilic agricultural biogas plant near Cologne, Germany. The samples were collected in a monthly interval a month over a period of 587 days and analysed based on 16S amplicon sequences by Martin Fisher et al.. We used three exemplanary samples from this dataset and accomplished and analysed the process from raw reads to MAGs. 
Link of the dataset:
https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/1751-7915.13313

## From raw reads to MAGs 

For the assembly and analyse of Metagenome Assembled Genomes (MAGs) we will follow these steps:

1. Pre-processing the raw reads (trimming adapters, removing chimeras, removing phiX174 virus sequences…)
2. Assemble reads into contigs/fasta
3.  Asses quality of assemblies
4.  Bin contigs into MAGs
5.  Asses completeness, contamination, and strain heterogeneity of your MAGs

## Tools used in the course: 
|Tool |	Version |	Repository|
----------------------------
|fastqc |	0.11.9 |	FastQC|
|fastp|	0.22.0	|fastp|
|megahit	|1.2.9|	megahit|
|samtools|	1.9|	samtools|
|QUAST|	5.0.2|	quast|
|Bowtie2|	2.4.5|	bowtie2|
|Concoct|	1.1.0|	CONCOCT|
|MetaBAT2|	2.12.1|	Metabat2|
|DASTool|	1.1.5|	DAS_Tool|
|anvi´o	|7.1	|anvi’o|
|GUNC|	1.0.5|	GUNC|

## HPC system and batch script

### HPC
Most of the steps were run on a High Performance Computing-Server (HPC) via Bash scripts. 
Therefore, we needed to access the HPC-Server of the CAU:
```
ssh -X sunamXXX@caucluster.rz.uni-kiel.de
```
### Batch script

The scripts were written in Sublime Text and executed though the Linux Shell via
```
sbatch name.sh
```
The batch script should contain the following information:

|Parameter|	
-------------
|#SBATCH |	Slurm batch script directive|
|--partition= or -p	|Slurm partition (~batch class)|
|--job-name= or -J	|Job name|
|--output= or -o	|Stdout file|
|--error= or -e	|Stderr file; if not specified, stderr is redirected to stdout file|
|--nodes= or -N |	Number of nodes|
|--ntasks-per-node=	|Number of tasks per node; number of MPI processes per node|
|--cpus-per-task= or -c|	Number of cores per task or process|
|--mem=<size[units]>|	Real memory required per node; default unit is megabytes (M); use G for gigabytes|
|--time= or -t|	Walltime in the format "hours:minutes:seconds"|


All packages and programs needed are already installed into one conda environment.
To activate the environment:
```
conda activate /home/sunam226/.conda/envs/anvio
```
or 
```
source activate /home/sunam226/.conda/envs/anvio
```

Example bash script:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=example
#SBATCH --output=example.out
#SBATCH --error=example.err
#SBATCH --partition=all
#SBATCH --reservation=biol217
#activate environment
module load miniconda3/4.7.12.1
conda activate /home/sunam226/.conda/env/anvio
#Commandlines
example.command
#this prints the required resources into your logfile
jobinfo
```

# 1. Pre-processing the raw reads

Creating our working directory
```
mkdir /work_beegfs/sunam228/Day2
```
Copying the raw reads into our working directory
```
cp /home/sunam226/Day2/0_raw_reads/*.fastq.gz /work_beegfs/sunam228/day2
cd ./day2
```
Quality control of the raw reads using **FASTQC**

```
for i in *.gz; do fastqc $i -o output_folder/; done
```

Depending on the phred-score, the sequences were trimmed with **Fastp**
```
for i in `ls *_R1.fastq.gz`;
do
 second=`echo ${i} | sed 's/_R1/_R2/g'`
 fastp -i ${i} -I ${second} -R _report -o ../clean_reads/"${i}" -O
../clean_reads/"${second}" -t 6 -q 20
done
```

# 2. Assembling the reads into contigs

The processed (clean) reads were assembled into contigs using **Megahit**
```
cd /work_beegfs/sunam228/Day2/clean_reads
megahit -1 BGR_130305_R1.fastq.gz -1 BGR_130527_R1.fastq.gz -1
BGR_130708_R1.fastq.gz -2 BGR_130305_R2.fastq.gz -2 BGR_130527_R1.fastq.gz -2
BGR_130708_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 -o
../megahit -t 20
```

For visualization, the contigs were transformed into FASTG files and opened in **Bandage**
```
megahit_toolkit contig2fastg 99 final.contigs.fa > final.contigs.fastg
```
### `Question 1`

Submit your generated figure and explain shortly what you can see

![image](https://user-images.githubusercontent.com/123310950/218768317-2cf3f1c7-1e33-4871-a9b8-5dba3fb5e419.png)

The picture shows the de novo created contigs from the processed/clean reads

- There are no bubbles, branches or unresolved ends in the assembly, but many short contigs
- The colour is selected based on assumed MAG assignment (bins)

# 3.  Asses quality of assemblies

The quality assessment of the assemblies from **megahit** was performed with **Quast**. The output is given as PDF and
html.
```
cd /work_beegfs/sunam228/3_coassembly/
metaquast -t 6 -o ../3_metaquast/ -m 1000 final.contigs.fa
```
### `Question 2`

- What is your N50 value?
  - 2963

- Why is this value relevant?
  - The N50 defines assembly quality in terms of contiguity

- How many contigs are assembled?
  - 57414
- What is the total length of the contigs?
  - 145675865

# 4.  Bin contigs into MAGs

Binning of the contigs:
```
anvi-script-reformat-fasta final.contigs.fa -o
/work_beegfs/sunam228/Day3/contigs.anvio.fa --min-len 1000 --simplify-names--
report-file name_conversion.txt
```
Afterwards the bins were mapped in a loop, so every final.contigs.fasta recieved a corresponding mapping file
```
module load bowtie2
cd ./2_fastp/
for i in `ls *mapped_R1.fastq.gz`;
do
 second=`echo ${i} | sed 's/_R1/_R2/g'`
 bowtie2 --very-fast -x ../4_mapping/contigs.anvio.fa.index -1 ${i} -2
${second} -S ../4_mapping/"$i".sam
done
```
The output are sequence mapping files (**.sam**), which were converted to binary alignment and map files
(**.bam**)
```
module load samtools

for i in *.sam; do samtools view -bS $i > "$i".bam; done
```
## 4.1 Anvi'o
The Ananlysis and Visualization plattform for microbial 'Omics combines many of the computational
strategies of data-enabled microbiology.

**convqqase** = anvi’o contigs-db database that contains key information associated with your sequences

```
anvi-gen-contigs-database -f ./4_mapping/contigs.anvio.fa -o ./5_anvio_profiles/contigs.db -n 'biol217'
```

The command will:
- Compute k-mer frequencies for each contig
- Soft-split contigs longer than 20,000 bp into smaller ones
- Identify open reading frames using Prodigal, the bacterial and archaeal gene finding program

Next HMM search on the contigs

HMM = "Basically, in anvi’o, `H`idden `M`arkov `M`odels (or `HMMs` for short) are used to search for specific genes (like SCGs) with known functions in a larger dataset"
```
anvi-run-hmms -c ./5_anvio_profiles/contigs.db
```
## 4.2 Access anvi’o interactive

To access (everytime) :
```
 srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=all /bin/bash
 ```
 ```
source activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

anvi-display-contigs-stats contigs.db
```

Open new Terminal
```
ssh -L 8060:localhost:8080 sunam###@caucluster.rz.uni-kiel.de
ssh -L 8080:localhost:8080 node###
```
Open http://127.0.0.1:8060/ in your browser

This program shows you simple stats of your contigs database that may help you not only assess your assembly output, but also estimate the number of bacterial and archaeal genomes to recover.

![image](https://user-images.githubusercontent.com/123310950/218780754-232ab5b6-1b77-46f8-9d4f-39eec1f8d9cb.png)
![image](https://user-images.githubusercontent.com/123310950/218780803-db21d8dc-45f0-4523-bc6b-3fd8c3851a5d.png)

## 4.3 Binning with Anvi'o

Sorted and indexed the .bam files with **samtools** in anvi'o.
```
for i in *.bam; do anvi-init-bam $i -o ../5_anvio_profiles/"$i".sorted.bam; done
```

### Creating an Anvi’o profile

An anvi’o profile stores sample-specific information about contigs. Profiling a BAM file with anvi’o using anvi-profile creates a single profile that reports properties for each contig in a single sample based on mapping results.
```
mkdir ./profiling/
cd ./5_anvio_profiles/
for i in `ls *.sorted.bam | cut -d "." -f 1`; do anvi-profile -i "$i".bam.sorted.bam.sorted.bam -c ./contigs.db -o ../profiling/$i; done
```
This command line will return a folder where you will find the following files:

- RUNNLOG.txt the anvi´o log of your run
- PROFILE.db An anvi’o database that contains key information about the mapping of short reads from multiple samples to your contigs

The created folder containins profiles.db and a .txt log

Processing of contigs will include:

- The recovery of mean coverage, standard deviation of coverage, and the average coverage for the inner quartiles (Q1 and Q3) for a given contig
- The characterization of single-nucleotide variants (SNVs) for every nucleotide position

### Merge the profiles coming from your different samples into one profile:
```
anvi-merge ./6_profiling/BGR_130305/PROFILE.db ./6_profiling/BGR_130527/PROFILE.db ./6_profiling/BGR_130708/PROFILE.db -o ./6_profiling/merged_profiles -c ./5_anvio_profiles/contigs.db --enforce-hierarchical-clustering
```
When you runanvi-merge:

It will merge everything and create a merged profile
It will attempt to create multiple clusterings of your splits using the default clustering configurations.

### Genome Binning 

Basically, in metagenomic binning, you’re trying to group together a bunch of contigs that all belong to the same genome using various metrics like tetranucleotide frequency, differential coverage, completion, etc .“

In this course two different Binners were used:
- Metabat2 
- Concoct

The result of both Binners were consolidated using **DASTool**

#### Binning with Metabat2
```
anvi-merge ./BGR_130305_mapped_R1/PROFILE.db ./BGR_130527_mapped_R1/PROFILE.db
./BGR_130708_mapped_R1/PROFILE.db -o ./ -c ../5_anvio_profiles/contigs.db --
enforce-hierarchical-clustering
```
#### Binning with CONCOCT
```
anvi-cluster-contigs -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db
-C consolidated_bins --driver dastool -T 20 --search-engine diamond -S
METABAT,CONCOCT --log-file log_consolidation_of_bins --just-do-it
anvi-summarize -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db -o
/PATH/TO/SUMMARY_consolidated_bins -C consolidated_bins
```

#### Consolidating bins with DASTool
```
anvi-cluster-contigs -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db
-C consolidated_bins --driver dastool -T 20 --search-engine diamond -S
METABAT,CONCOCT --log-file log_consolidation_of_bins --just-do-it
anvi-summarize -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db -o
/PATH/TO/SUMMARY_consolidated_bins -C consolidated_bins
```

### `Question 3`

Number of Archaea bins you got from MetaBAT2?
- 3

Number of Archaea bins you got from CONCOCT?
- 2

Number of Archaea bins you got after consolidating the bins?
- 2


### MAGs Quality Estimation

Visualizing and evaluating the results.
Estimate your genomes completeness and contamination levels. You can assess the quality of your bins by using

```
anvi-estimate-genome-completeness -c ./contigs.db -p ./5_anvio_profiles/merged_profiles/PROFILE.db -C consolidated_bins > genome_completeness_das
tool-txt
```
To check what collections you generated you can use:
```
anvi-estimate-genome-completeness -p ./5_anvio_profiles/merged_profiles/PROFILE.db -c ./contigs.db --list-collections
```

COLLECTIONS FOUND

METABAT (48 bins, representing 8939 items).
CONCOCT (156 bins, representing 57626 items).
consolidated_bins (41 bins, representing 13825 items).

#### Open anvi'o interactive
```
srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=all /bin/bash
```
```
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

anvi-interactive -p ./5_anvio_profiles/merged_profiles/PROFILE.db -c ./contigs.db -C METABAT
```
Open New Terminal
```
ssh -L 8060:localhost:8080 sunam228@caucluster-old.rz.uni-kiel.de
ssh -L 8080:localhost:8080 node077
```
Open in a Browser: http://127.0.0.1:8060/

Anvi-interactive gives you the possibility to manually inspect and work on bins.
![image](https://user-images.githubusercontent.com/123310950/218786622-cd94617e-0d70-407d-a701-4794b5110ab4.png)

### `Question 4`:

Which binning strategy gives the best quality archea bins:
- dastool (consolidated) gives the best quality, which was expected since it it not a binner itself and just
merges the results of the binners. 
Metabat2 seems to be better than CONCOCT for this data.

How many archaea bins do you get of high quality
- 2

### Bin refinement

From this point only the archea bins were used. 

A summary folder was created with anvi-summarize:
```
anvi-summarize -c ../contigs.db -p ./merged_profiles/PROFILE.db -C consolidated_bins -o ./summary --just-do-it
```

The summarize folder containins a comprehensive overview of the collection and statistics created by anvio.

Then, the archea bins were copied to a seperate folder

```
cd ./summary/bin_by_bin

mkdir ../../ARCHAEA_BIN_REFINEMENT

cp ./summary/bin_by_bin/Bin_Bin_1_sub/*.fa ./ARCHAEA_BIN_REFINEMENT/

cp ./summary/bin_by_bin/Bin_METABAT__25/*.fa ./ARCHAEA_BIN_REFINEMENT/

```
# 6.  Asses completeness, contamination, and strain heterogeneity of your MAGs


