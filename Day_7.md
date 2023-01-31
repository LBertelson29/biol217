# Day 7

# RNA-Seq analysis (Transcriptomics) - Tutorial

## Aim

The aim of this tuorial is to train you in:

- Setting up RNA-Seq data analysis pipeline
- RNA-Seq data pre-processing
- RNA-Seq Data Analysis
- Differential Gene Expression
- Data Visualization
- Functional enrichment Analysis

You will:

- Download the required files
- Setting up the files in specific locations
- Checking read quality
- Align the reads to a reference sequence
- Calculate the coverage
- Perform gene wise quantification
- Calculate differential gene expression
- Much more....

## Tools used

Tool 	 	

- READemption 2.0.3 	
- DESeq2 4.2 	
- edgeR 3.32 	
- limma 3.46 
- R 4.1.0 	
- RStudio 1.4.1717 
- Kallisto 0.48.0 


## Dataset to be used in the example

The dataset you will use today comes from a publication by Prasse et al. 2017.
How to download the data to be used for RNA_seq Analysis?

1. Activate the environment:
```
conda activate /home/sunam226/.conda/envs/grabseq
```
2. Download the data specifying the SRA:
3. 
`SRA Number: 	SRP081251`

`Open the paper from this Prasse et al. 2017, find out the SRR numbers, quantity of samples and treatments, and write down here`:

[link](https://www.tandfonline.com/doi/full/10.1080/15476286.2017.1306170)

Accessible through GEO Series accession number GSE85456 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85456). A shell script that can be used to reproduce the RNA-seq analysis can be retrieved from Zenodo at https://zenodo.org/record/59989 (DOI: 10.5281/zenodo.59989).

Search for NCBI or accesssion. It will show the link for the National Library of Medicine/

There you can find the SRR number, # of Spots, # of Bases, Size and the Date when they were published. 

```
grabseqs -t 4 -m metadata.csv SRR***
```

`SRR`: 
1. SRR4018516
2. SRR4018515
3. SRR4018514
4. SRR4018517

### Nevigate to new folder: 
```
bash mkdir fastq cd fastq 
``` 
### Download the data specifying the SRA: 

```
grabseqs -t 4 -m metadata.csv SRR4018514 SRR4018515 SRR4018516 SRR4018517
```
## Batch-Script for grabseqs
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=grabseq
#SBATCH --output=grapseq.out
#SBATCH --error=grapseq.out
#SBATCH --partition=all
#SBATCH --reservation=biol217

#load your anvio environment (path needs to be adjusted)

module load miniconda3/4.7.12.1

source activate /home/sunam226/.conda/envs/grabseq

#navigate to working directory
cd /work_beegfs/sunam228/Day7/fastq

grabseqs -t -m ./metadata.csv SRR4018514 SRR4018515 SRR4018516 SRR4018517
```

`--> Does not work because of grabseqs`

## Alternative: 

```
fasterq-dump SRR 
```

## Batch-script for fasterq-dump
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=fasterq-dump
#SBATCH --output=fasterq-dump.out
#SBATCH --error=fasterq-dump.out
#SBATCH --partition=all
#SBATCH --reservation=biol217

#load your anvio environment (path needs to be adjusted)

module load miniconda3/4.7.12.1

source activate /home/sunam226/.conda/envs/grabseq

#navigate to working directory
cd /work_beegfs/sunam228/Day7/fastq

fasterq-dump SRR 4018514 
fasterq-dump SRR 4018515 
fasterq-dump SRR 4018516
fasterq-dump SRR 4018517
```
`Script does not work`

### Do it manual 

```
fasterq-dump SRR4018514 SRR4018515 SRR4018516 SRR4018517
```


Try the dry run for grabseqs

```
grabseqs sra -l -t 4 -m metadata.csv SRP081251

```

`Rename the fastq-files`

sRNA_R1 or 2 for the mutant

wt_R1 or 2 for the wildtype


## You can also download it manually from the NIH
fasta and fastq-files, filtered or clipped


## Download the SRP from Kröger et al. 2013
```
grabseqs sra -l -t 4 -m metadata.csv SRP013870
```


# Quality Control: How to run fastqc?

1. Activate the environment:
   
```
module load fastqc
fastqc -t 4 -o fastqc_output *.fastq.gz
```

## Batch-script for fastqc 
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.out
#SBATCH --partition=all
#SBATCH --reservation=biol217

#load your anvio environment (path needs to be adjusted)

module load miniconda3/4.7.12.1
source activate /home/sunam226/.conda/envs/grabseq

#navigate to working directory
cd /work_beegfs/sunam228/Day7/fastq/

fastqc -t 4 -o fastqc_output *.fastq
```


## Compine all folders and sub-folders into one with this command: 

```
multiqc -d . -o multiqc_output 
```

For my direction:
```
multiqc -d . -o ./multiqc_output/ 
```
-d will also count sub folders

2. Open the fastqc_output folder and check the quality of the reads.

via html file


## Notes from lecture
### `Important files for RNA-seq analysis: FASTA, FASTQ and GFF/GTF`

## For a experiment: 
Deeper sequencing or more replicates?--> more replicates

How many replicates? --> Depends on the money, but 3 are good

Read length? --> 50bp for single and 100bp for paired, not over 300bp.

Why to we use paired-ends? because we know the distance (between 15bp or 34bp)

# READemption

- To evaluate this dataset you will use READemption "a pipeline for the computational evaluation of RNA-Seq data."
- We have already installed the pipeline in your CAUCLUSTER IDs, however, you can read this documentation for more details on installation procedure.
- As the analysis will take a while, you will run a script that includes a pipeline with all READemption commands needed for this analysis.


## The dataset you will use today comes from a publication by Kröger et al. 2013.

The aim of the study was to "present a simplified approach for global promoter identification in bacteria using RNA-seq-based transcriptomic analyses of 22 distinct infection-relevant environmental conditions. Individual RNA samples were combined to identify most of the 3,838 Salmonella enterica serovar Typhimurium promoters" (Kröger et al. 2013). Twenty two 22 different environmental conditions were used to study the effect on gene expression.

Here you will use a subset of the original dataset including two replicates from two conditions:

1. InSPI2 "an acidic phosphate-limiting minimal media that induces Salmonella pathogenicity island (SPI) 2 transcription".
This is suspected to create an environmental shock that could induce the upregulation of specific gene sets.
   
2. LSP growth in Lennox Broth medium

## Workflow

## The READemption pipeline is divided into 5 steps:

1. create - Create a project folder and the required subfolders

`example:rna_seq`

create environment
```
conda activate /home/sunam226/.conda/envs/reademption

```
look at the softwares in this environment: 
```
conda list
```

Create folders:
```
reademption create --project_path READemption_analysis --species salmonella="Salmonella Typhimurium"
```

It creates 2 Folders in READemption_analysis
Input and Output
In Input we habe 3 folders
There we will add our 3 files: FASTA, FASTQ and GFF/GFT
FASTA in reads, annotations and reference



1. Prepare the reference sequences and annotations

Save `Ftp source` one time, to be used later several times:
```
FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645/
```

Download the .fasta files for genome:
```
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa $FTP_SOURCE/NC_016810.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa $FTP_SOURCE/NC_017718.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa $FTP_SOURCE/NC_017719.fna
wget -O READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa $FTP_SOURCE/NC_017720.fna
```

## Now, we will download the genome annotation file:
```
wget -P READemption_analysis/input/salmonella_annotations https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.gff.gz
```

unzip the file:
```
gunzip READemption_analysis/input/salmonella_annotations/GCF_000210855.2_ASM21085v2_genomic.gff.gz
```
`Manually: Remove in the fasta files in the first line <gi|386730646|ref|>`

or with: 

Change the headers of fasta files .fa with the header of annotation file, using following commands:
```
sed -i "s/>/>NC_016810.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_016810.fa
sed -i "s/>/>NC_017718.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017718.fa
sed -i "s/>/>NC_017719.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017719.fa
sed -i "s/>/>NC_017720.1 /" READemption_analysis/input/salmonella_reference_sequences/NC_017720.fa
```
Now we will download the raw reads as mentioned here:
```
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/InSPI2_R2.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R1.fa.bz2
wget -P READemption_analysis/input/reads http://reademptiondata.imib-zinf.net/LSP_R2.fa.bz2
```

## After downloading the files make a .bash script and submit the job, as follows:

## Batch-script READemption pipeline
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


# Our Project: Methanosarcina mazei


1. create - Create a project folder and the required subfolders

`project_publications`

1. create environment
```
conda activate /home/sunam226/.conda/envs/reademption

```

1. Create project folder:
```
reademption create --project_path READemption_analysis --species Mmazei="Methanosarcina mazei"
```

3. Open NCBI and search Methanosarcina mazei

Download the FASTA and GFF-files


4. (Change the headers of fasta files .fa with the header of annotation file, using following commands) 
   `already changed`

## Sbatch-script READemption pipeline
```
#!/bin/bash
#SBATCH --job-name=reademption_M
#SBATCH --output=reademption_M.out
#SBATCH --error=reademption_M.err
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
reademption align --fastq -p 4 --poly_a_clipping --project_path READemption_analysis
reademption coverage -p 4 --project_path READemption_analysis
reademption gene_quanti -p 4 --features CDS,tRNA,rRNA --project_path READemption_analysis
reademption deseq -l sRNA_R1.fastq,sRNA_R2.fastq,wt_R1.fastq,wt_R2.fastq -c sRNA,sRNA,wt,wt -r 1,2,1,2 --libs_by_species Mmarzei=sRNA_R1,sRNA_R2,wt_R1,wt_R2 --project_path READemption_analysis
reademption viz_align --project_path READemption_analysis
reademption viz_gene_quanti --project_path READemption_analysis
reademption viz_deseq --project_path READemption_analysis
conda deactivate
jobinfo
```

