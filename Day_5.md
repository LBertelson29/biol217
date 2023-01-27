# Day 5
# MAGs Quality Estimation
Estimate your genomes completeness and contamination levels.
You can assess the quality of your bins by using
```
anvi-estimate-genome-completeness -c /PATH/TO/contigs.db -p merged_profiles/PROFILE.db -C consolidated_bins
```
`For my path`

```
anvi-estimate-genome-completeness -c ./contigs.db -p ./5_anvio_profiles/merged_profiles/PROFILE.db -C consolidated_bins
```

### Visualizing and evaluating the results 

To check what collections you generated you can use:

```
 anvi-estimate-genome-completeness -p merged_profiles/PROFILE.db -c /PATH/TO/contigs.db --list-collections
 ```
`For my path` 
```
anvi-estimate-genome-completeness -p ./5_anvio_profiles/merged_profiles/PROFILE.db -c ./contigs.db --list-collections
```

 ### `Open anvi'o interactive`

```
srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=all /bin/bash
```
- node077 -
```
conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

anvi-interactive -p ./5_anvio_profiles/merged_profiles/PROFILE.db -c ./contigs.db -C METABAT
```
`Open New Terminal`
```
ssh -L 8060:localhost:8080 sunam228@caucluster-old.rz.uni-kiel.de
```
```
ssh -L 8080:localhost:8080 node077
```

`Open google chrome and paste`

http://127.0.0.1:8060 or http://127.0.0.1:8080


Anvi-interactive gives you the possibility to manually inspect and work on bins.
- you can set all parameters that you want 

![image](Pictures/Collection__consolidated_bins__for_merged_profiles.svg)




## `Questions` 



Which binning strategy gives you the best quality for the Archaea bins?

How many Archaea bins do you get that are of High Quality? 



