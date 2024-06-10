#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=140G
#SBATCH -p hpcbioamd
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jholmes5@illinois.edu
#SBATCH -J pileup
#SBATCH --array=1-2
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/src/slurm-out

#Author - Jessica Holmes (github @jrkirk61), HPCBio, Carver Biotechnoloy Center, University of I
llinois

# The -D option above it to put the slurm out files in their own directory.

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load popscle/20210504-IGB-gcc-8.2.0

### load sample names
basename=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" ../Ey_Op_basenames.txt)

### specifying the output folder
cd ../../results/

### Run popscle freemuxlet
popscle freemuxlet --plp freemuxlet/${basename}/ \
  --out freemuxlet/${basename}/freemuxed \
  --nsample 3

