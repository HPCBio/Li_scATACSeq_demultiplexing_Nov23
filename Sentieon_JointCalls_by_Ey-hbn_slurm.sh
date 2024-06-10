#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p hpcbio
#SBATCH --mem 200G
#SBATCH --mail-type=END,FAIL
#SBATCH -J Sentieon_JointCalls_by_group
#SBATCH --mail-user=jholmes5@illinois.edu
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/src/slurm-out/
#SBATCH -o JointCalls_by_group.ou.%a
#SBATCH -e JointCalls_by_group.er.%a

#Author - Gloria Rendon (github @grendon), HPCBio, Carver Biotechnoloy Center, University of Illinois

# Only the variables INPUTDIR, OUTPUTDIR, should need to be changed from sample to sample. Below that, the "sentieon driver" code should also be edited to merge the exact files you want to merge. 

set -x 
echo `date`
echo `hostname`

set -euo pipefail


module load sentieon/202308


echo `date`
echo Declare variables

INPUTDIR=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/results/VariantCalling/Ey_1/Variants/
OUTPUTDIR=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/results/VariantCalling/Ey_1/JointCalls/
GENOME=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/data/references/dm6_genome.fna
DICT=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/data/references/dm6_genome.dict
INDEX=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/data/references/dm6_genome.fna.fai

echo `date`
echo Sanity check

if [ ! -s $GENOME ]
then
echo $GENOME file not found
exit 1
fi
if [ ! -d $INPUTDIR ]
then
echo $INPUTDIR path not found
exit 1
fi

if [ ! -d $OUTPUTDIR ]
then
mkdir -p $OUTPUTDIR
fi


cd $OUTPUTDIR
ln -s $INPUTDIR/* ./

echo `date`
echo start Joint Calls for Ey_hbn

sentieon driver -r ${GENOME} --algo GVCFtyper \
-v Ey_hbn_1_S3_L006_VARIANT.GVCF  -v Ey_hbn_1_S3_L007_VARIANT.GVCF  \
-v Ey_hbn_2_S4_L006_VARIANT.GVCF -v Ey_hbn_2_S4_L007_VARIANT.GVCF \
Samples_Pool1_Run_JOINT_RAW_VARIANTS.VCF 

if [ ! -s Samples_Pool1_Run_JOINT_RAW_VARIANTS.VCF ]
then

echo GVCFtyper failed Ey_hbn
exit 1
fi


echo `date`
echo Done Joint Calls

