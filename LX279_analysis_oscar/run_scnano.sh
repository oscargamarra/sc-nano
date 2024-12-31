#! /bin/bash
#SBATCH --parsable -p cpu --job-name=elnano --time=148:00:00 --gres=tmpspace:50G --mem=40GB --nodes=1 --mincpus=8
conda activate /hpc/pmc_holstege/oscar/anaconda3/envs/scNanoGPS
P_DIR="/hpc/pmc_holstege/oscar/scNanoGPS/snvcalling"
FASTQ="/hpc/pmc_holstege/scg/nanopore/pod5/LX279/LX279-sup.fastq.gz"
REF_GENOME="/hpc/pmc_holstege/rstudio/oscar_nanopore/snv_calling/genome.fa"
IND_GENOME="/hpc/pmc_holstege/rstudio/oscar_nanopore/snv_calling/genome.fa.mmi"
ncores=8
ANNOVAR="/hpc/pmc_holstege/oscar/annovar/"
ANNOVAR_DB="/hpc/pmc_holstege/oscar/annovar/humandb/"
ANNOVAR_GV="hg38"
ANNOVAR_XREF="/hpc/pmc_holstege/oscar/annovar/humandb/hg38_refGene.txt"

python3 $P_DIR/../scanner.py -i $FASTQ -t $ncores --scanning_region=130
python3 $P_DIR/../assigner.py -t $ncores
python3 $P_DIR/../curator.py -t $ncores --ref_genome $REF_GENOME --idx_genome $IND_GENOME
python3 $P_DIR/../reporter_SNV.py -t $ncores --ref_genome $REF_GENOME --annovar $ANNOVAR --annovar_db $ANNOVAR_DB --annovar_gv $ANNOVAR_GV
