# Day 6

# Pangenomics - comparing genomes



In this tutorial we will combine both the previously assembled MAGs and reference genomes for a phylogenetic and functional genome comparison. This tutorial follows the workflow of the anvi'o miniworkshop and the pangenomics workflow.

Workflow:

- Recap on the Batch Script
- Evaluating the contigs databases
- Create pangenome from individual bins/genomes
- Compare the data phylogenetically (ANI)
- Visualizing the pangenome
- Interpreting and ordering the pangenome
- BONUS: BlastKoala

Folder Structure: 02_contigs-dbs, 03_pangenome

Program or Database 	      Function

anvi'o = Wrapper for genome comparissons

DIAMOND = creates high-throughput protein alignments

pyANI = calculates genome similarities based on average nucleotide identity

BlastKOALA = Onlinetool which creates metabolic networks for a given genome, based on the KEGG database

KEGG = Kyoto Encyclopaedia of Genes and Genomes (Database)

NCBI COG = Clusters of Orthologous Genes (Database)


## 1. A recap on the batch script and for loops

To create a batch script copy the dummy from here, or one of your older scripts

cp ....sh .

The batch script should contain:

1. The shebang
2. Processing requirements as #SBATCH commands
    - Reservation
    - Nodes to use
    - CPUs (for multithreading)
    - memory requirements
    - time
    - working directory
    - log files
    - partitions
3. Stdin and Stderr paths
4. Slurm modules needed for the task
5. The command
6. jobinfo

```
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=1:30:00

#SBATCH -D ./
#SBATCH --output=./3_fasta-for-anvio.out
#SBATCH --error=./3_fasta-for-anvio.out
#SBATCH --partition=all

# for pangenome
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

# set working directory by navigating there
cd ....

# Insert your command here


# provides information on resource requirements as stdout
jobinfo
```

# Start of Pangenomics

## 2. Evaluation of our starting databases 
### (Directory: 02_contigs-dbs)


New set of contigs.dbs
They contain:
- MAGs from the Biogasreactor
-  and a complete Methanogen Genome

In order to get us started we will visualize and compare these bins in a summary overview.


This is done with the function:
```
#get direct access to a HPC compute node
srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 /bin/bash

#activate the conda environment
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

# start anvi'o interactive display
anvi-display-contigs-stats *db
```
`Open new Terminal`
```
ssh -L 8060:localhost:8080 sunam228@caucluster-old.rz.uni-kiel.de
ssh -L 8080:localhost:8080 node010
```
`Open a chrome browser` and enter the following IP. http://127.0.0.1:8060

Task: Take some time to click through the views and compare the MAGs. 
Add a screenshot of your output to your documentation. 

![Image](Pictures/Contigs_DB_Stats_1.png)
![Image](Pictures/Contigs_DB_Stats2.png)
![Image](Pictures/Contigs_DB_Stats3.png)

## `Question`
### How do the MAGs compare in size and number of contigs to the full genome?

The number and size of MAGs differ in compare to the full genome trough all bins. 
The number of contigs in Methano_Mflavescens is 1, as for the bins the number differ from 164 to 334. 
The total length for Methano_Mflavescens is 3,283,688, in the bins the total length is between 1,336,430 and 2,637,263.


### Based on the contig numbers, sizes and number of marker genes (HMM hits), which two MAGs are the best:

based on contig number : Methano_Bin13 & Methano_Bin8 (lowest contig numbers)

based on size : Methano_Bin13 & Methano_Bin9 (highest size)

based on number of marker genes : Methano_Bin13 & Methano_Bin3 / Methano_Bin8 (most identical HMM hits)

Methano_bin13 and Metano_bin09 (highest size)
### and which is the worst?
based on contig number : Methano_Bin5 (highest contig number)

based on size : Methano_Bin10 (smallest size)

based on number of marker genes : Methano_Bin5 least identical HMM hits)




When done, close the window and ```Ctrl+C``` in the command lines and ```exit```

# 3. Making a Pangenome 
## (Directory: 03_pangenome)

A pangenome visualizes entire genomes for comparisson.

It can show essential and accessory gene clusters, phylogenetic relationships and genome qualities.

## 3.1 Create an external genomes file

To tell the programm which genomes and MAGs it should use, we will create the "external genomes file".

The external genomes file contains one column with the genome/bin name and one with its path (where it is saved).

We already have a folder with all the genome databases we want to compare (02_contigs-dbs). Anvi'o has a script to create the input information for us:

TASK: Complete the following line, and use it on the login node.

```
anvi-script-gen-genomes-file --input-dir ? -o external-genomes.txt
```

My line: 
```
anvi-script-gen-genomes-file --input-dir ../02_contigs-dbs -o external-genomes.txt
```

Now look into your file to verify whether it looks accurate.

`Tip` use cat or head
 - cat: First 10 lines
 - tail: last 10 lines

## `Output`

# 3.2 Estimate genome completeness

To avoid any nasty suprises by adding a bad bin or incomplete genome to the pangenome, estimate genome completeness. This will give you information on the quality of your MAGs and genomes.

Question: The command provides its output as a table to the standard output of the terminal. What can you add to the code to direct output to, e.g. a .txt file?

```
anvi-estimate-genome-completeness -e external-genomes.txt > genome-completeness.txt
```
We want to specifically look at redundancy and completeness.

## `Question`
How do the bins compare to isolate genomes? Would you remove one, based on the output of the completeness estimation?

# 3.3 Remove unwanted genomes (Directory: 02_contigs-dbs)

As we have some MAGs with a low completion, we will remove them from our pangenome. Common practice is, to consider only genomes with > 70% completion and < 10% redundancy.

For this go back to 02_contig-dbs, create a new directory "discarded" and mv the "bad MAGs_dbs" to this folder.

```
mkdir discarded
```
```
mv /work_beegfs/sunam228/Day6/02_contigs-dbs/Bin10.db discarded
mv /work_beegfs/sunam228/Day6/02_contigs-dbs/Bin5.db discarded
```
Return to 03_pangenome and recreate the external genomes file.

```
anvi-script-gen-genomes-file --input-dir ? -o external-genomes-final.txt
```
`My line`
```
anvi-script-gen-genomes-file --input-dir ../02_contigs-dbs/ -o external-genomes-final.txt
```
# 3.4 Creating the pangenome database (Directory: 03_pangenome)

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

We will combine the next two steps in one BATCH script with the following computing requirements:

`change SBATCH Settings`: --nodes=1, --cpus-per-task=10, --mem=500M, --time=00:05:00

Look for the following commands and settings and complete this in your batch script:
 ```
 anvi-gen-genomes-storage -e ? -o ?

anvi-pan-genome -g ? --project-name ? --num-threads 10
```
My Batch-script: 
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


# 4. Genome similarity based on average nucleotide identity (ANI) 
## (Directory: 03_pangenome)

The next step calculates the genome similarity to each other. The most commonly used approach is average nucleotide identity using the MUMmer algorithm to align each genome. The result of this is used as a measure to determine how related the genomes are and whether you have discovered a new species. Usually the cutoff for the species boundary is set at 95-96% identity over a 90% genome coverage [Ciufo, et al., 2018; Jain, et al. (2018)].

Once anvi'o has calculated the genome similarity, you can use its output to organize your genomes based on their relatedness.

Depending on the amount of genomes you are using, this step can be quite memory intensive.

Find out what the following parameters mean and complete the command in a BATCH script:

`SBATCH Settings`: --nodes=1, --cpus-per-task=10, --mem=600M, --time=00:02:00
```
anvi-compute-genome-similarity --external-genomes ? --program ? --output-dir ? --num-threads ? --pan-db ?
```
`My SBATCH-script`
```

```
