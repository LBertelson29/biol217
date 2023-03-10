# biol217 (Marine)Microbial Omics - from sample to function

Aim of the course:
- Practical and theoretical knowledge of a bioinformatic workflow
- Assembly and analyse of MAGs 
- Skills in the command line interface
- documentation of bioinformatic analyses

# 1. Metagenome Assembled Genomes (MAGs)

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
|------|--------|------------|
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

|Parameter|	  |
|-------------|-------------|
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

## 1.1 Pre-processing the raw reads

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

## 1.2 Assembling the reads into contigs

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

**Submit your generated figure and explain shortly what you can see**

![image](https://user-images.githubusercontent.com/123310950/218768317-2cf3f1c7-1e33-4871-a9b8-5dba3fb5e419.png)

The picture shows the de novo created contigs from the processed/clean reads

- There are no bubbles, branches or unresolved ends in the assembly, but many short contigs
- The colour is selected based on assumed MAG assignment (bins)

## 1.3  Asses quality of assemblies

The quality assessment of the assemblies from **megahit** was performed with **Quast**. The output is given as PDF and
html.
```
cd /work_beegfs/sunam228/3_coassembly/
metaquast -t 6 -o ../3_metaquast/ -m 1000 final.contigs.fa
```
### `Question 2`

**What is your N50 value?**
  - 2963

**Why is this value relevant?**
  - The N50 defines assembly quality in terms of contiguity

**How many contigs are assembled?**
  - 57414
  
**What is the total length of the contigs?**
  - 145675865

## 1.4  Bin contigs into MAGs

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
## 1.4.1 Anvi'o
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
## 1.4.2 Access anvi’o interactive

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

## 1.4.3 Binning with Anvi'o

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

**Number of Archaea bins you got from MetaBAT2?**
- 3

**Number of Archaea bins you got from CONCOCT?**
- 2

**Number of Archaea bins you got after consolidating the bins?**
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

**Which binning strategy gives the best quality archea bins:**
- dastool (consolidated) gives the best quality, which was expected since it it not a binner itself and just
merges the results of the binners. 
Metabat2 seems to be better than CONCOCT for this data.

**How many archaea bins do you get of high quality?**
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
## 1.5 Asses completeness, contamination, and strain heterogeneity of your MAGs
## GUNC
Use GUNC to check run chimera detection (Chimeric genomes are genomes wrongly assembled out of two or more genomes coming from separate organisms).

Genome UNClutter (GUNC) is “a tool for detection of chimerism and contamination in prokaryotic genomes resulting from mis-binning of genomic contigs from unrelated lineages.”

to use GUNC , activate the following environment:
```
conda activate /home/sunam226/.conda/envs/gunc
```
Process all your files in one run:

```
cd ../../archea_bin_refinement
mkdir GUNC
cd /work_beegfs/sunam228/Day5/5_anvio_profiles/archea_bin_refinement
for i in *.fa; do gunc run -i "$i" -r /home/sunam226/Databases/
gunc_db_progenomes2.1.dmnd --out_dir GUNC --threads 10 --detailed_output; done
```

To see if a MAG is chimeric or not (clade seperation score) the output file (.tsv) was opened.

### `Questions 5`
**Do you get bins that are chimeric?**
- The Metabat-bins are non-chimeric
- The Concoct-bins from kingdom to class are chimeric

**In your own words (2 sentences max), explain what is a chimeric bin:**
- chimeric contigs: Wrongly assembled contigs from genomic fragments.
- chimeric bins: contigous fragments from different sources are sorted in the same genomic bin.

## Manual bin refinement

As large metagenome assemblies can result in hundreds of bins, pre-select the better ones for manual refinement, e.g. > 70% completeness.

Before you start, make a copy/backup of your unrefined bins the ARCHAEA_BIN_REFINEMENT (manual refinement will overwrite the unrefined files).

Activate Anvi'o interactive

```
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1
```
```
srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --nodelist=node002 /bin/bash
node010
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1
```
```
anvi-refine -c ../contigs.sb -C consolidated_bins -p ./merged_profiles/PROFILE.db --bin-id Bin_METABAT__25
```
Open new terminal
```
ssh -L 8060:localhost:8080 sunam228@caucluster-old.rz.uni-kiel.de
ssh -L 8080:localhost:8080 node010
```
Open in a browser: http://127.0.0.1:8060/

Sort your bins by GC content, by coverage or both

The interface allows you to categorize contigs into separate bins (selection tool)
- Unhighlighted contigs are removed when the data is saved
- Evaluate taxonomy and duplicate single copy core genes
- Remove contigs

For refinement use clustering based on only differential coverage, and then only based on sequence composition in search for outliers

### `Questions 6`

**Does the quality of your Archaea improve? (hint: look at completeness, strain heterogeneity)**

- The completeness of the bins decreases (from 97.4 to 93.4)
- The strain heterogeneity stays the same

Table + output figure
![image](https://user-images.githubusercontent.com/123310950/218795797-8f590230-d8cb-474c-b930-d635a726b9bd.png)

### `Question 7`

**How abundant are the archea bins in the 3 samples**
- Metabat: 1.76 | 1.14 | 0.58
- Concoct: 0.96 | 0.00 | 0.40

## Taxonomic assignment
Add taxonomic annotations to your MAG:

```
anvi-run-scg-taxonomy -c ../contigs.db -T 20 -P 2
```
**Anvi-estimate-scg-taxonomy** -This program makes quick taxonomy estimates for genomes, metagenomes, or bins stored in your contigs-db using single-copy core genes (SCGs)

```
anvi-estimate-scg-taxonomy -c ../contigs.db --metagenome-mode
```
The abundance of rRNAs in the dataset was estimated:
```
anvi-estimate-scg-taxonomy -c ./4_mapping/contigs.db -p
./5_anvio_profiles/merged_profiles/PROFILE.db --metagenome-mode --compute-scgcoverages --update-profile-db-with-taxonomy
```
The output was saved from the terminal:
```
anvi-estimate-scg-taxonomy -c ./4_mapping/contigs.db -p
./5_anvio_profiles/merged_profiles/PROFILE.db --metagenome-mode --compute-scgcoverages --update-profile-db-with-taxonomy > temp.txt
```
Final summary to get comprehensive info about your consolidated bins:
```
anvi-summarize -p /PATH/TO/merged_profiles/PROFILE.db -c /PATH/TO/contigs.db -o /PATH/TO/SUMMARY_consolidated_bins -C consolidated_bins
```

### `Questions 8`

**Did you get a species assignment to the Archaea bins previously identified?**

**Concuct:**

total_length = 1948187
num_contigs = 1054
N50 = 1883
GC_content = 43.0927763589817
percent_completion = 73.6842105263158
percent_redundancy = 1.31578947368421
t_domain = Archaea
t_phylum = Halobacteriota
t_class = Methanosarcinia
t_order = Methanosarcinales
t_family = Methanosarcinaceae
t_genus = Methanosarcina
t_species = Methanosarcina flavescens


**Metabat:**

total_length = 1859269
num_contigs = 250
N50 = 8819
GC_content = 59.5651181815842
percent_completion = 97.3684210526316
percent_redundancy = 5.26315789473684
t_domain = Archaea
t_phylum = Halobacteriota
t_class = Methanomicrobia
t_order = Methanomicrobiales
t_family = Methanoculleaceae
t_genus = Methanoculleus
t_species = Methanoculleus sp012797575


Does the HIGH-QUALITY assignment of the bin need revision?
- Yes.



# 2. Pangenomics - comparing genomes
## 2.1 Evaluation of the starting databases

In this tutorial we will combine both the previously assembled MAGs  from the biogas reactor and reference genomes for a phylogenetic and functional genome comparison.

**Tools:**

- anvi'o = Wrapper for genome comparissons
- DIAMOND = creates high-throughput protein alignments
- pyANI = calculates genome similarities based on average nucleotide identity
- BlastKOALA = Onlinetool which creates metabolic networks for a given genome, based on the KEGG database
- KEGG = Kyoto Encyclopaedia of Genes and Genomes (Database)
- NCBI COG = Clusters of Orthologous Genes (Database)

To compare the bins a summary overview was vizualized.
New set of contigs.dbs
They contain:
- MAGs from the Biogasreactor
- Complete Methanogen Genome

In order to get us started we will visualize and compare these bins in a summary overview:
```
srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 /bin/bash
```
```
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1
```
```
anvi-display-contigs-stats *db
```
Open new Terminal
```
ssh -L 8060:localhost:8080 sunam228@caucluster-old.rz.uni-kiel.de
ssh -L 8080:localhost:8080 node010
```

Open in a Browser: http://127.0.0.1:8060/

![Image](Pictures/Contigs_DB_Stats_1.png)
![Image](Pictures/Contigs_DB_Stats2.png)
![Image](Pictures/Contigs_DB_Stats3.png)

### `Question` 9
**How do the MAGs compare in size and number of contigs to the full genome?**

-The size of MAGs differ in compare to the full genome trough all bins. The size of the MAGs are smaller than the full Methano_mflavenscens genome.
-The number of contigs is much higher in the MAG bins (164 to 334), since the full genome only cosists of 1 contig.
-The longest contigs vary between 14 kb and 65 kb, compared to 3.28 mb of the full genome of Methano_mflavenscens.
-Single copy core genes (HMM) differ from 43 to 80. the complete genome contains 72. 

**Based on the contig numbers, sizes and number of marker genes (HMM hits), which two MAGs are the best:**

- based on contig number : Methano_Bin13 & Methano_Bin8 (lowest contig numbers)
- based on size : Methano_Bin13 & Methano_Bin9 (highest size)
- based on number of marker genes : Methano_Bin13 & Methano_Bin3 / Methano_Bin8 (most identical HMM hits)

Methano_bin13 and Metano_bin09 (highest size)

**And which is the worst?**
- based on contig number : Methano_Bin5 (highest contig number)
- based on size : Methano_Bin10 (smallest size)
- based on number of marker genes : Methano_Bin5 least identical HMM hits)


## 2.2 Making a Pangenome 

A pangenome visualizes entire genomes for comparisson.

It can show essential and accessory gene clusters, phylogenetic relationships and genome qualities.

### Creating external genome files

To tell the programm which genomes and MAGs it should use, we will create the "external genomes file".

The external genomes file contains one column with the genome/bin name and one with its path (where it is saved).

 Anvi'o has a script to create the input information for us:
 ```
 anvi-script-gen-genomes-file --input-dir ../02_contigs-dbs -o external-genomes.txt
 ```
Now look into your file to verify whether it looks accurate:
```
cat external_genomes.txt 
name contigs_db_path
Methano_Bin1 /work_beegfs/sunam228/Day6/02_contigs-dbs/Bin1.db Methano_Bin10
/work_beegfs/sunam228/Day6/02_contigs-dbs/Bin10.db Methano_Bin13
/work_beegfs/sunam228/Day6/02_contigs-dbs/Bin13.db Methano_Bin3
/work_beegfs/sunam228/Day6/02_contigs-dbs/Bin3.db Methano_Bin5
/work_beegfs/sunam228/Day6/02_contigs-dbs/Bin5.db Methano_Bin8
/work_beegfs/sunam228/Day6/02_contigs-dbs/Bin8.db Methano_Bin9
/work_beegfs/sunam228/Day6/02_contigs-dbs/Bin9.db Methano_Mflavescens
/work_beegfs/sunam218/Day6/02_contigs-dbs/Mflavescens.db
```
### Estimate genome completeness

To receive information about the quality of the MAGs, the genome completeness of each MAG was estimated trough the factors redundancy and completeness. 
```
anvi-estimate-genome-completeness -e external_genomes.txt -o ./genomecompleteness.txt
```

### `Question 10`

**How do the bins compare to isolate genomes? Would you remove one, based on the output of the completeness estimation?**
- Some bins deviate  greatly to the isolate genome. The bins show a lower completion and a higher redundancy compared to isolated genomes.
- I would remove Methano_Bin5 and Methano_Bin10 based on the completeness estimation. 

## Remove unwanted genomes

As we have some MAGs with a low completion, we will remove them from our pangenome. Common practice is, to consider only genomes with **> 70% completion** and **< 10% redundancy**.

The unwanted MAGs (> 70% completion and < 10% redundancy) were moved to a new folder "discarded". 

```
mkdir discarded
```
```
mv /work_beegfs/sunam228/Day6/02_contigs-dbs/Bin10.db discarded
mv /work_beegfs/sunam228/Day6/02_contigs-dbs/Bin5.db discarded
```

A new external-genomes output was created for the updated bins:
```
anvi-script-gen-genomes-file --input-dir ../02_contigs-dbs/ -o external-genomes-final.txt
```
## 2.3 Creating the pangenome databas

In anvi'o we will need to generate two artifacts (similar to when working with assemblies)

The first is the **genomes-storage.db**, which corresponds to an individual contigs.db, but merges all individual genomes you are working with into one database. 

The database contains:

1. all genome fasta files
2. the gene annotations (HMMs, SCGs) which were added before
3. any new annotations and genome comparisons we will make

The second file is the **pan-genome.db**. It is similar to the profile you generate to annotate your bins.

This will contain:

1. genome similarities based on gene amino acid sequences.
2. resolved gene clusters
3. any post-analysis of gene clusters, downstream analyses and visualisations

The next two steps were combined in one BATCH script:

```
#!/bin/bash

#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=500M
#SBATCH --time=00:05:00
#SBATCH --job-name=pangenome_database
#SBATCH -D ./
#SBATCH --output=./pangenome_database.out
#SBATCH --error=./pangenome_database.out
#SBATCH --partition=all

# for pangenome
source activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

# set working directory by navigating there
cd /work_beegfs/sunam228/Day6/03_pangenome

anvi-gen-genomes-storage -e external-genomes-final.txt -o ./Methano_GENOMES.db

anvi-pan-genome -g Methano_GENOMES.db --project-name Methano_pangenome --num-threads 10
```

## 2.4 Genome similarity based on average nucleotide identity (ANI) 

In anvi'o we will need to generate two artifacts, similar to when working with assemblies. The first is the genomes-storage.db, which corresponds to an individual contigs.db, but merges all individual genomes you are working with into one database. The files themselves will be a bit leaner, than all files together, making it easier to share and publish those.

The database contains:

1. all genome fasta files
2. the gene annotations (HMMs, SCGs) which were added before
3. any new annotations and genome comparisons we will make

The second file is the pan-genome.db. It is similar to the profile you generate to annotate your bins.

This will contain:

1. genome similarities based on gene amino acid sequences.
2. resolved gene clusters
3. any post-analysis of gene clusters, downstream analyses and visualisations

We will combine the next two steps in one Batch script with the following computing requirements:

```
#!/bin/bash

#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=500M
#SBATCH --time=00:05:00
#SBATCH --job-name=pangenome_database
#SBATCH -D ./
#SBATCH --output=./pangenome_database.out
#SBATCH --error=./pangenome_database.out
#SBATCH --partition=all

# for pangenome
source activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

# set working directory by navigating there
cd /work_beegfs/sunam228/Day6/03_pangenome

anvi-gen-genomes-storage -e external-genomes-final.txt -o ./Methano_GENOMES.db

anvi-pan-genome -g Methano_GENOMES.db --project-name Methano_pangenome --num-threads 10
```


# 2.4 Genome similarity based on average nucleotide identity (ANI) 

The next step calculates the genome similarity to each other. The most commonly used approach is average nucleotide identity using the MUMmer algorithm to align each genome. The result of this is used as a measure to determine how related the genomes are and whether you have discovered a new species. Usually the cutoff for the species boundary is set at **95-96% identity over a 90% genome coverage**.


```
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=600M
#SBATCH --time=00:02:00
#SBATCH --job-name=ANI
#SBATCH -D ./
#SBATCH --output=./ANI.out
#SBATCH --error=./ANI.out
#SBATCH --partition=all

# for pangenome
source activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

# set working directory by navigating there
cd /work_beegfs/sunam228/Day6/03_pangenome

anvi-compute-genome-similarity --external-genomes external-genomes-final.txt --program pyANI --output-dir ANI  --num-threads 10 --pan-db ./Methano_pangenome/Methano_pangenome-PAN.db
```
# 2.5 Visualizing the pangenome 

 Access to a HPC compute node:

```
srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --reservation=biol217 --partition=all /bin/bash

#activate the conda environment
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

anvi-display-pan -p Methano_pangenome-PAN.db -g ../Methano_GENOMES.db -P 8083
```

### `Questions 11`

**Which INPUT FILES do we need?**

- Input files from the pangenome analysis.

  -p PAN_DB, --pan-db PAN_DB
                        Anvi'o pan database (default: None)
  -g GENOMES_STORAGE, --genomes-storage GENOMES_STORAGE
                        Anvi'o genomes storage file (default: None)



**Task**
**Write the command and use the additional flag -P.**
```
anvi-display-pan -p Methano_pangenome-PAN.db -g ../Methano_GENOMES.db -P 8083
```

**What is the -P flag for?**

-P INT, --port-number INT
                        Port number to use for anvi'o services. If nothing is
                        declared, anvi'o will try to find a suitable port
                        number, starting from the default port number, 8080.
                        (default: None)

Open new terminal in MobaXterm and start the tunnel: 
```
ssh -L 8060:localhost:8083 sunam228@caucluster-old.rz.uni-kiel.de 
```
```
ssh -L 8083:localhost:8083 node010
```
Open browser: http://127.0.0.1:8060

# 2.6 Interpreting and ordering the pangenome (interactive interface)
**TASKS**: Genome similarity

1. Remove combined homogeneity, functional homogeneity, geometric homogeneity, max num parsimonay, number of genes in gene cluster and number of genomes gene cluster has hits from the active view. `Tip`: Play with Height

2. Create a "Bin-highlight" including alls SCGs and name it accordingly. 

3. Cluster the genomes based on Frequency

### `Questions 12`
**Based on the frequency clustering of genes, do you think all genomes are related? Why?**
- Based on the SCGs, which are clustered in the same region and only occur once in every genome, all genomes should classify as Archea.
Bin9 seems to be closer related to Mflavescens. The two share more gene clusters than the other bins, which are more similar to each other that to Bin09 and M.Flavescens. 

**How does the reference genome compare to its closest bin?**
- The reference genom and the genom from Bin09 are pretty simular, with mostly identical gene clustering, especially in the SCG region. 

**What ranges are used determine a prokaryotic species?**
- Prokaryotic species cut off 95%. 

So Bin9 and refence seem to be one species. The other four differ from the reference, but share a species with each other
Differences within the clustered bins appear at 99%
Higher percentage=higher approximate relatedness. 95% cutoff for the same species in prokaryotes.

**How high can you go until you see changes in ANI?**
- Changes in AI appear at 99%.
With 99.5% cutoff, Bin01 is not related to the other Bins.
99.6% cutoff = Bin13 related to 3 and 8, Bin8 and 3 loosely related
99.7% cutoff= only relatedness between Bin9 and the reference.

**What does the ANI clustering tell you about genome relatedness?**
- A higher percentage indicates a higher approximate relatedness. 

**TASKS**: Functional Profiling

- Using the Search Function, highlight all genes in the KEGG Module for Methanogenesis
- Create a new bin called "Methanogenesis" and store your search results in this bin.


### `Question 13`: 
**How are Methanogenesis genes distributed across the genome?**
- They are distributed across the entire genome, no clusters recognizable.

**Task**
- Google COG Categories and select one you are interesed in. Create a new bin, find your Category in the Pangenome and add it to this selection.


 FINAL VIEW 

![Image](Pictures/Methano_pangenome(1).svg)


### `Question 14`: 
**What observations can you make regarding the geometric homogeneity between all genomes and the functional homogeneity?**
- 20 Genes match among the bins and are located in the marked SCGs region. 
The AA structure differs though the genes. 
The most similarities in the AA make-up can be found between the reference and Bin09, with simililar SNPs and deletions/insertions.


 ![Image](Pictures/Screenshot%202023-01-30%20at%2015-20-40%20GC_00000248%20detailed.png)


### BONUS: BlastKoala

Outside of anvi'o there are a range of tools available to investigate your organisms metabolism. One of these is BlastKOALA, which generates a metabolic profile of your genome based on the KEGG database.

**Task**: Check out the BlastKOALA Results for this Methanogen.

### `Question 15`
**Can the organism do methanogenesis? Does it have genes similar to a bacterial secretion system?**
- Yes, the organism can do methanogenesis. 7 similar genes were found, which are used in the bacterial secretion system.

# 3. RNA-Seq analysis (Transcriptomics)

The aim: 

- Setting up RNA-Seq data analysis pipeline
- RNA-Seq data pre-processing
- RNA-Seq Data Analysis
- Differential Gene Expression
- Data Visualization
- Functional enrichment Analysis

**Tools** 	 	

- READemption 2.0.3 	
- DESeq2 4.2 	
- edgeR 3.32 	
- limma 3.46 
- R 4.1.0 	
- RStudio 1.4.1717 
- Kallisto 0.48.0 

**Dataset**
The dataset you will use today comes from a publication by Prasse et al. 2017.
https://www.tandfonline.com/doi/full/10.1080/15476286.2017.1306170

## 3.1 Downloading the data to be used for RNA_seq Analysis:

To download the data used in the paper the SRR number was used:
```
conda activate /home/sunam226/.conda/envs/grabseq
grabseqs -t 4 -m ./metadata.csv SRR4018514 SRR4018515 SRR4018516 SRR4018517
```
Download the entire project:
```
conda activate /home/sunam226/.conda/envs/grabseq
grabseqs -t 4 -m SRP081251
```

Or download the data manually (fasta and fastq-files, filtered or clipped).

### `Question 15`
**How to find the SRR number from a paper?**
Search in the paper directly for SRR, NCBI or Accession with ctrl+f
Open the link
Look for a SRA or project numer
Click on each samples and look for SRR-Number 


## 3.2  Quality control (fastqc & fastp)

1. Run fastqc
```
module load fastqc
fastqc -t 4 -o fastqc_output *.fastq.gz
```
2. A fastqc output folder was created
3. Run fastp
4. Multiqc was used to see a summary of all files at once
```
mkdir ../qc_reports
for i in *.fastq.gz; do fastqc -t 4 -o ../qc_reports/fastqc_output $i; done
for i in *.fastq.gz; do fastp -i $i -o ${i}_cleaned.fastq.gz -h
../qc_reports/${i}_fastp.html -j ${i}_fastp.json -w 4 -q 20 -z 4; done
cd ..
multiqc -d . -o ./qc_reports/multiqc_output 
```

## 3.3 READemption

- **READemption** = A pipeline for the computational evaluation of RNA-Seq data
- Pipeline in the CAUCLUSTER IDs

The READemption pipeline is divided into 5 steps:

**1. Creating a project folder and the required subfolders**
```
conda activate /home/sunam226/.conda/envs/reademption
reademption create --project_path READemption_analysis --species
salmonella="Salmonella Typhimurium"
```

**This will create 2 Folders in READemption_analysis**
- Input and Output

**2. Prepare the reference sequences and annotations**

Save `Ftp source for easier access:

FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645/

2.2 Download the .fasta files for genome:

```
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa $FTP_SOURCE/NC_016810.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa $FTP_SOURCE/NC_017718.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa $FTP_SOURCE/NC_017719.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa $FTP_SOURCE/NC_017720.fna
```

2.3 Download the genome annotation file (.gff3) and unzip it
```
wget -P READemption_analysis/input/salmonella_annotations https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.gff.gz

gunzip READemption_analysis/input/salmonella_annotations/GCF_000210855.2_ASM21085v2_genomic.gff.gz
```

2.4 Change the headers 

Change the headers of the fasta files to the header of annotation file (annotation file + reference file
need to look the same: Sequence ID first):

By hand: Remove everything infront of the NZ-number in sublime text editor, than safe.

Via command:
```
sed -i "s/>/>NC_016810.1 /"
READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa
sed -i "s/>/>NC_017718.1 /"
READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa
sed -i "s/>/>NC_017719.1 /"
READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa
sed -i "s/>/>NC_017720.1 /"
READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa
```

2.5 Downloading the raw reads

```
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R2.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R2.fa.bz2
```

## 3.4 Run READemption pipeline

```
#!/bin/bash
#SBATCH --job-name=reademption
#SBATCH --output=reademption.out
#SBATCH --error=reademption.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --qos=long
#SBATCH --time=0-05:00:00
#SBATCH --partition=all
#SBATCH --export=NONE
#SBATCH --reservation=biol217


#activating conda
module load miniconda3/4.7.12.1
source activate /home/sunam226/.conda/envs/reademption
reademption align -p 4 --poly_a_clipping --project_path READemption_analysis
reademption coverage -p 4 --project_path READemption_analysis
reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis
reademption deseq -l InSPI2_R1.fa.bz2,InSPI2_R2.fa.bz2,LSP_R1.fa.bz2,LSP_R2.fa.bz2 -c InSPI2,InSPI2,LSP,LSP -r 1,2,1,2 --libs_by_species salmonella=InSPI2_R1,InSPI2_R2,LSP_R1,LSP_R2 --project_path READemption_analysis
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis
conda deactivate
jobinfo
```
## 3.5 IGB
### 3.5.1 Working with the Integrated Genome Browser (IGB)
There you can visualize:
bam files
fasta-files
gff
wig 

Load for sequence: fasta-file
Load for annotation: gff-file
Load for coverage : bw-files

**Bacterial Ribo-seq data analysis by the HRIBO pipeline**
![Image](Pictures/Bacterial%20Ribo-seq%20data%20analysis%20by%20the%20HRIBO%20pipeline.%20(A)%20In.png)


### 3.5.2 Searching genes in IGB
**Task**: Find the ORF rncoding csrA and some of it features

Location of the gene csrA : NC_016856.1:2,991,176-2,991,380

![Image](Pictures/IGB_csrA.png)


The locus tag for the gene encoding CsrS : STM14_3412

1. Start codon (canonical ATG or alternative GTG,TTG) : ATG
2. Its stop codon (TAG,TGA,TAA) : TAA
3. Its length in amino acids: 61 (-stop codon) 
4. Its SD (consensus AGGAGG at -7 to -4): AGGAG at -7
5. The name of the upstream(before) gene: STM14_3413/alaS
6. Do you think csrA is translated? Why?: Yes, as the coverage of the Ribo-seq covers the whole gene (the coverage fits the gene). UTRs clearly visible.

With the riborex-file you can see if the gene is translated or not and more information. 

**Task**: Find a significantly differential expressed gene betweeen wt and scrA-mutant 

## `Gene gltK`

![Image](Pictures/Gene_gltK.png)


## 4. R - Introduction

# R-introduction

1. Start with basic codes

```
getwd()  #shows the working directory (path)

setwd() #set working directory

dir.create("data") #creates folder with the name data

dir.create("data/raw_data") # creates folder (raw_data) inside a folder (data)

# Learning basic codes ----------- 

# to set a variable 

x <- 2+2  (the same like this x=2+2, but better)


# what type is my variable? 
x>-5
class(x)

-->"numeric"

x <- "Hello World"
class (x)

--> "charachter"

x<- TRUE
class(x)

--> "logical"

x <- 1+2i
class(x)

--> "colplex"

x <- charToRaw ("Hello World")
class(x)

-->"raw"
x <-data.frame(matrix(1:6,nrow = 2, ncol = 3))
class(x)

--> "data.frame"


6=6 # 6 is the variable 6

6==6 #6 is 6

6!=6 # 6 is not 6 

my_obj <-40 

data() # shows data already in R

View(data) # shows the data

class(dataname$variable) # shows what kind of variable this is

class(iris$Species) 

"factor"
```
2. Creating plots 


```
plot(iris)

boxplot(data=iris, x=iris$Petal.Length)


boxplot(data=iris, iris$Petal.Length~iris$Species)

stripchart()

hist()

qqnorm()

pairs()
```
3. Installing packages for better plots (for example ggplot2)
```
go to packages --> install --> choose the one you want from CRAN
```
or 
```
install.packages("ggplot2")

install.packages("tidyverse")

```
install 2 packages (or more) at once
```
install.packages(c("readxl", "plotly"))

```

4. Ggplot 2 

Help function 
```
?ggplot # Help function
```

For plotting, we need: 
- mapping
- aesthetics
- geom

```
ggplot(data = iris, mapping = aes(x = iris$Species, y = iris$Sepal.Length)) + geom_boxplot()
```
or 

```
ggplot(data = iris, mapping = aes(x = Species, y = Sepal.Length)) + geom_boxplot()

```
change color

fill = color in aes
```
ggplot(data = iris, mapping = aes(x = Species, y = Sepal.Length, fill = Species)) + geom_boxplot()
```
or with col = 
```
ggplot(data = iris, mapping = aes( Sepal.Length, Petal.Length, col = Species)) + geom_point()
```
Change the shape of the datapoints with shape = 
```
ggplot(data = iris, mapping = aes( Sepal.Length, Petal.Length, shape = Species)) + geom_point()
```

5. Saving a plot 
First assign a plot name
```
plot1 <-ggplot(data = iris, mapping = aes( Sepal.Length, Petal.Length, size = Petal.Length)) + geom_point()
```
save it 

as PDF
```
plot1 + ggsave("plot1.pdf", height = 6, width = 8, units = "in", dpi = 300)
```
or tiff
```
plot1 + ggsave("plot1.tiff", height = 6, width = 8, units = "in", dpi = 300)

```
or png 
```
plot1 + ggsave("plot1.png", height = 6, width = 8, units = "in", dpi = 300)

```
`dpi = dots per inch`
(like pixels)

compress an image (for mailing etc., only for png or tiff)

```
plot1 + ggsave("plot1.tiff", height = 6, width = 8, units = "in", dpi = 300. compression = "lzw")
```


### `Question 16`

How to reshape the Data?

from long to wide 
```
tidyr::gather(df1, key= 'class values', value = 'numeric values')
```


6. Normal-distribution 

With histogram, only for numeric values
```
hist(iris$Sepal.Length)
```
7. Converting table into data.frame

trees<-data.frame(tress)

8. Learning plotting with own choosen dataset
```
Trees<-data.frame(trees)

tree1<-ggplot(data = trees, mapping = aes(Height, Volume, col = Height)) + geom_boxplot()
tree2<-ggplot(data = trees, mapping = aes(Girth, Volume)) + geom_point()
tree3<-ggplot(data = trees, mapping = aes(Volume, Height, size = Girth, col = Girth)) + geom_point()
tree4<-ggplot(data = trees, mapping = aes(Girth, Volume)) + geom_line()
tree5<-ggplot(data = trees, mapping = aes(Girth, Height)) + geom_line()

```

![Image](Pictures/R/trees1.png)

![Image](Pictures/R/trees2.png)

![Image](Pictures/R/tree3.png)

![Image](Pictures/R/tree4.png)

![Image](Pictures/R/tree5.png)

8. Set the repo
```
setRepositories() # in console
```

9. Import dataset 

In Environment `Import Dataset`

10. What kind of plot can i use for my data? 

`Google andrew abela chart guide`

![Image](Pictures/chart_choose.png)



11. Plotting a heatmap 

11.1 Heatmap with name of the gene, sequence type and log2fold

```
library(readxl)
df <- read_excel("csrA-WT_sorted.xlsx", sheet = "Sheet6")
View(df)

library(ggplot2)

ggplot(df, aes(df$Name, df$seq_type, fill=df$log2fold_change)) + geom_tile()

ggplot(df, aes(df$seq_type, df$Name,fill=df$log2fold_change)) + geom_tile()
```
![Image](Pictures/Heatmap1.png)

**change color** (+ black outline)
```
ggplot(df, aes(df$seq_type, df$Name,fill=df$log2fold_change)) + geom_tile()+ scale_fill_gradient(low = "white", high = "red")+theme_linedraw() 
```
![Image](Pictures/Heatmap2.png)

**Changing the x,y and legend-titels** 

```
ggplot(df, aes(df$seq_type, df$Name,fill=df$log2fold_change)) + geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+theme_linedraw()+
  labs(x = "Sequence type", y = "Genname")+
  guides(fill=guide_legend(title="log2fold"))
```

![Image](Pictures/Heatmap3.png)


11.2 Heatmap with identifier, sequence type and lod2fold

first, gather the data 

```
library(tidyr)

df1<-tidyr::gather(df1, key= 'Identifier', value = 'log2FC')
```

plot the heatmap 
```
ggplot(df1, aes(df1$method, df1$Identifier,fill=df1$log2FC)) + geom_tile()+
  scale_fill_gradient(low = "white", high = "red")+theme_linedraw()+
labs(x = "Sequence type", y = "Identifier")+
  guides(fill=guide_legend(title="log2FC"))
```

![Image](Pictures/Heatmap4.png)







