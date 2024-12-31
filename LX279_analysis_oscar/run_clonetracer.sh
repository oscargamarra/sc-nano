#!/bin/bash
#SBATCH --parsable
#SBATCH -p gpu
#SBATCH --job-name=elnano
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:1,tmpspace:50G
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

python run_clonetracer.py -i /hpc/pmc_holstege/rstudio/oscar_clonetracer/pat0.json -n pat0-2 -o /hpc/pmc_holstege/rstudio/oscar_clonetracer -t 300 -g

