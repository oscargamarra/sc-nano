#! /bin/bash
#SBATCH --parsable -p cpu --job-name=elnano --time=24:00:00  --gres=tmpspace:50G --mem=2GB --nodes=1 --mincpus=1

#### RASUSA (V-0.6.1) is a tool for downsampling long read data
#For doc: https://github.com/mbhall88/rasusa
#reads function subsamples fastq files
#-s argument selects a seed for reproducibility. In this case 33
#-n argument selects the number of reads of downsampled fastq
#-b argument selects the number of base pairs of downsampled fastq
#-V gives version, for easier reproducibility

#First attempt downsampled to number of reads.
#Wasn't used, just left it here for the record
rasusa reads /hpc/pmc_holstege/scg/nanopore/pod5/LX279/LX279-sup.fastq.gz -s 33 -n 60471134 -o /hpc/pmc_holstege/rstudio/oscar_clonetracer/LX279_example.fastq

#Second attempt, downsampled to number of base pairs
#Used this one for remapping with epi2me
rasusa reads /hpc/pmc_holstege/scg/nanopore/pod5/LX279/LX279-sup.fastq.gz -o /hpc/pmc_holstege/rstudio/oscar_clonetracer/LX279-sup-downsampled.fastq -b 48681908640 -s 123 -V
