# Day 2

# Batch-script

```
!/bin/bash
SBATCH --nodes=1
SBATCH --cpus-per-task=4
SBATCH --mem=10G
SBATCH --time=1:00:00
SBATCH --job-name=fastqc
SBATCH --output=/work_beegfs/sunam226/fastqc.out
SBATCH --error=/work_beegfs/sunam226/fastqc.err
SBATCH --partition=all
SBATCH --reservation=biol217

#only change name in line 6-8

module load miniconda3/4.7.12.1

conda activate anvio
```

fastqc

#load your anvio environment (path needs to be adjusted)
```
source activate /home/sunamXXX/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1
```


#set environmental variables here (anything you want to run multiple times)
outdir=/work_beegfs/sunamXXX/out_dir
meta=/work_beegfs/sunamXXX/.txt
runname=Run2023-01-20
trunclength=240

#commands are printed into the log
set -xu

#navigate to working directory
cd $HOME

# write your command here, if you want to use some environemtal variables, these can be accessed via $nameofvariable

# HERE YOU WRITE YOUR COMMANDS
#
#
rm -rf /home/sunam226/.conda/envs/viromics

# Finish each script with this (prints Done and a Unicorn into your logile), then you know everything has run through
echo "Done removing"

printf '\U1F984\n'

#this prints the required resources into your logfile
jobinfo

urs@Kurs015:~/Desktop/biol217$ cd
kurs@Kurs015:~$ ssh -X sunam228@caucluster.rz.uni-kiel.de
sunam228@caucluster.rz.uni-kiel.de's password: 
Permission denied, please try again.
sunam228@caucluster.rz.uni-kiel.de's password: 
Last failed login: Tue Jan 24 13:46:08 CET 2023 from ukifam-c168054.ifam.uni-kiel.de on ssh:notty
There was 1 failed login attempt since the last successful login.
Last login: Mon Jan 23 11:34:02 2023 from ukifam-c168054.ifam.uni-kiel.de


[sunam228@caucluster2 ~]$ cd /home/sunam226 
[sunam228@caucluster2 sunam226]$ ls
anviscript  Databases  Day1  Day10  Day2  Day3  Day4  Day5  Day6  Day7  Day8  Day9  roadmap_bio217.png  test_day1
[sunam228@caucluster2 sunam226]$ cp anviscript /home/sunam228
[sunam228@caucluster2 sunam226]$ cd /home
[sunam228@caucluster2 home]$ ls

so komme ich in andere sunam um sachen zu downloaden 

squeue -u sunam 228 

da kann man sehen, was gerade passiert 

$ conda activate /home/sunam226/.conda/envs/anvio
falls was schief l??uft, das hier nehmen

so conda activate environment
conda activate /home/sunam226/.conda/envs/anvio

module load fastqc 

und dann $for i in *.gz; do fastqc $i -o output_folder/; done


$ squeue -u sunam228

$ fastp --h --> zeigt, was man alles schneiden kann von den strings 


sbatch ist ein Protokoll, was ich an sich ganz dann abschicken kann 

ls -l zeigt an, was in meinem Ordner ist 

pwd zeigt meinen Pfad an 

sbatch assembly_script.sh  -> so lasse ich mein batch-Skript ausf??hren. Muss in den Ordner, wo mein Skript ist und das dann aufs??rhren 

Batch skript mit .sh abspeichern 

