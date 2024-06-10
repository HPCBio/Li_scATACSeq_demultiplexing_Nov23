# Li_scATACSeq_demultiplexing_Nov23

Repository of code used for demultiplexing Xin Li's scATACSeq data from November 2023.

## Order of scripts w/o filtering variants

_Recommended method_

1. Sentieon_Ey-1_varcalling_slurm.sh
2. Sentieon_Ey-2_varcalling_slurm.sh
3. Sentieon_JointCalls_by_Ey-hbn_slurm.sh
4. Demuxafy's script: sort_vcf_same_as_bam.sh
5. popscle_pileup_Ey-freemuxlet.sh
6. popscle_freemuxlet_Ey.sh 

## Order of scripts w/ filtering variants

1. Sentieon_Ey-1_varcalling_slurm.sh
2. Sentieon_Ey-2_varcalling_slurm.sh
3. Sentieon_JointCalls_by_Ey-hbn_slurm.sh
4. filtered-scripts/filter_joint_calls.sh
5. Demuxafy's script: sort_vcf_same_as_bam.sh
6. filtered-scripts/popscle_pileup_Ey-freemuxlet-filteredVCF.sh
7. filtered-scripts/popscle_freemuxlet_Ey_filteredVCF.sh
