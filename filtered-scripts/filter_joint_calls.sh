#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p hpcbioamd
#SBATCH --mem 15G
#SBATCH --mail-type=END,FAIL
#SBATCH -J vcf-filtering
#SBATCH --mail-user=jholmes5@illinois.edu
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/src/slurm-out/
#SBATCH --array=1

#Author - Jessica Holmes (github @jrkirk61), HPCBio, Carver Biotechnoloy Center, University of Illinois

# Change directories
cd /home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/results/VariantCalling/Ey_1/JointCalls/filtered/

Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted.VCF

# Load modules
module load GATK/4.4.0.0-Java-17.0.6
module load BCFtools/1.17-IGB-gcc-8.2.0


# 1. Separate out SNPs

gatk SelectVariants \
    -V ../Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted.VCF \
    -select-type SNP \
    -O Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly.vcf


# 2. GATK recommended filtering

gatk VariantFiltration \
    -V Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly.vcf \
    --cluster-size 3 \
    --cluster-window-size 10 \
    -filter "QD < 10.0" --filter-name "QD10" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_filtered.vcf


# 3. Remove filtered entries

grep "#" Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_filtered.vcf > \
	Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_filtered_PASSonly.vcf
grep -P "\tPASS\t" Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_filtered.vcf >> \
	Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_filtered_PASSonly.vcf


# 4. Remove non-exon variants
## Only needed for RNA-Seq data! 
 ##Compress VCF files w/ bgzipi
 #bgzip Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_filtered_PASSonly.vcf

 ##Index VCF file
 #bcftools index Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_filtered_PASSonly.vcf.gz

 ## Pull out exon variants
 #bcftools view -R ../../../../../data/references/clownfish_custom_exons.bed \
 #-o Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_PASSonly_exonsOnly.vcf \
 #Samples_Pool1_Run_JOINT_RAW_VARIANTS_bashsorted_SNPsOnly_filtered_PASSonly.vcf.gz

