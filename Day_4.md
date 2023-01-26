# Day_4

## Contigs data preparation

convqqase = anvi’o contigs-db database that contains key information associated with your sequences

```
anvi-gen-contigs-database -f contigs.anvio.fa -o contigs.db -n 'biol217'
```

## Batch-script for this part 
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=anvi_contigs
#SBATCH --output=anvi_contigs.out
#SBATCH --error=anvi_contigs.err
#SBATCH --partition=all
#SBATCH --reservation=biol217

#load your anvio environment (path needs to be adjusted)

module load miniconda3/4.7.12.1
source activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#navigate to working directory
cd /work_beegfs/sunam228/Day3/

anvi-gen-contigs-database -f ./4_mapping/contigs.anvio.fa -o ./5_anvio_profiles/contigs.db -n 'biol217'

``````
When you run this command, anvi-gen-contigs-database will 

   * Compute k-mer frequencies for each contig
  *  Soft-split contigs longer than 20,000 bp into smaller ones
   * Identify open reading frames using Prodigal, the bacterial and archaeal gene finding program



Next `HMM` search on the contigs

`HMM` = "Basically, in anvi’o, `H`idden `M`arkov `M`odels (or HMMs for short) are used to search for specific genes with known functions in a larger dataset"

```
anvi-run-hmms -c contigs.db
```
## Batch-script for HMM

```#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=hmm
#SBATCH --output=hmm.out
#SBATCH --error=hmm.err
#SBATCH --partition=all
#SBATCH --reservation=biol217

#load your anvio environment (path needs to be adjusted)

module load miniconda3/4.7.12.1
source activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#navigate to working directory
cd /work_beegfs/sunam228/Day3/

anvi-run-hmms -c ./5_anvio_profiles/contigs.db
```

Once you have your contigs database ready, and optionally your HMMs are run, you can take a quick look at it using the program `anvi-display-contigs-stats`


First you need to access anvi’o interactive 


To access (everytime) : 
```
 srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=all /bin/bash
```
  - `node068` - (remember!)
```
source activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

anvi-display-contigs-stats contigs.db
```

`Open new Terminal`

```
ssh -L 8060:localhost:8080 sunam228@caucluster.rz.uni-kiel.de
ssh -L 8080:localhost:8080 node###
```
`Open google chrome and paste`

http://0.0.0.0:8081

http://127.0.0.1:8060 or http://127.0.0.1:8080 (Remenber your server information)

This program shows you simple stats of your contigs database that may help you not only assess your assembly output, but also estimate the number of bacterial and archaeal genomes to recover.


## Questions Day_3 part 2

Number of bins you got from MetaBAT2?

Number of bins you got from CONCOCT?

Number of bins you got after consolidating the bins?

# Binning with ANVI´O

`ANVI´O`, an `AN`alysis and `VI`sualization platform for microbial ´`O`mics

