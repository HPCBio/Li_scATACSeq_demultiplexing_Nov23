#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p hpcbioamd
#SBATCH --mem 200G
#SBATCH --mail-type=END,FAIL
#SBATCH -J Li_Sentieon_DNA_varcalling
#SBATCH --mail-user=jholmes5@illinois.edu
#SBATCH -D /home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/src/slurm-out/
#SBATCH --array=3
#ignoreSBATCH -o DNA_varcalling_log_sample.ou.%a
#ignoreSBATCH -e DNA_varcalling_log_sample.er.%a

#Author - Gloria Rendon (github @grendon), HPCBio, Carver Biotechnoloy Center, University of Illinois

# Only the variables OUTPUTDIR, FILES, DBPATH, and SCRATCH should need to be changed from person to person, and/or sample to sample. Example of samplesheet.csv is also included in this repo.

set -x 
echo `date`
echo `hostname`

set -euo pipefail

THREADS=24
OUTPUTDIR=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/results/VariantCalling/Ey_1/
FILES=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/results/VariantCalling/ey_1_samplesheet.csv
SAMPLEID=`head -n $SLURM_ARRAY_TASK_ID $FILES | tail -n 1 | cut -d ',' -f1`
READ1=`head -n $SLURM_ARRAY_TASK_ID $FILES | tail -n 1 | cut -d ',' -f2`
READ2=`head -n $SLURM_ARRAY_TASK_ID $FILES | tail -n 1 | cut -d ',' -f3 | tr -d '\n' | tr -d '\r'`
GENOME=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/data/references/dm6_genome.fna
DICT=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/data/references/dm6_genome.dict
INDEX=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/data/references/dm6_genome.fna.fai
ADAPTERS=/home/groups/hpcbio/RNA-Seq/projects/lixin/2023Nov-scRNASeq/results/VariantCalling/adapter.fna
DBPATH=/scratch/jholmes5/Li_VarCalling
DBSNP=
SCRATCH=/scratch/jholmes5/Li_VarCalling/$SAMPLEID
SM=$SAMPLEID
ID=$SAMPLEID
PL=Illumina
#RLEN=150

echo `date`
echo Sanity check

if [ ! -s $GENOME ]
then
echo $GENOME file not found
exit 1
fi

if [ ! -s $READ1 ]
then
echo $READ1 file not found
exit 1
fi

if [ ! -s $READ2 ]
then
echo $READ2 file not found
exit 1
fi

if [ ! -d $OUTPUTDIR/trim ]
then
mkdir -p $OUTPUTDIR/trim
fi
if [ ! -d $OUTPUTDIR/sorted ]
then
mkdir -p $OUTPUTDIR/sorted
fi
if [ ! -d $OUTPUTDIR/dedup ]
then
mkdir -p $OUTPUTDIR/dedup 
fi
if [ ! -d $OUTPUTDIR/Variants ]
then
mkdir -p $OUTPUTDIR/Variants
fi
if [ ! -d $OUTPUTDIR/QC ]
then
mkdir -p $OUTPUTDIR/QC
fi
if [ ! -d $DBPATH ]
then
mkdir -p $DBPATH
fi
if [ ! -d $SCRATCH ]
then
mkdir -p $SCRATCH
fi


echo `date`
echo stage files to scratch $SCRATCH

SAMPLER1=$(basename $READ1 | tr -d '\n' )
SAMPLER2=$(basename $READ2 | tr -d '\n' )
REFERENCE=${DBPATH}/$(basename $GENOME | tr -d '\n' )



if [ ! -s $REFERENCE ]
then
cp ${GENOME}* $DBPATH/
cp $DICT $DBPATH/
cp $INDEX $DBPATH/
fi

cd $SCRATCH


if [ ! -s $SAMPLER1 ]
then
ln -s $READ1  $SAMPLER1
ln -s $READ2  $SAMPLER2
fi

echo `date`
echo start TRIMMING with sample ${SAMPLEID}
module load fastp/0.20.0-IGB-gcc-4.9.4

fastp -i ${SAMPLER1} \
	-I ${SAMPLER2} \
	-g \
	-w ${THREADS} \
	-q 10 \
	-o ${SAMPLEID}.R1.trimmed.fastq.gz \
	-O ${SAMPLEID}.R2.trimmed.fastq.gz \
	-h ${SAMPLEID}_fastp.html \
	-j ${SAMPLEID}_fastp.json \
	-l 25 \
  --adapter_fasta=$ADAPTERS 

module unload fastp/0.20.0-IGB-gcc-4.9.4
echo `date`
echo end TRIMMING with sample ${SAMPLEID}


if [ ! -s ${SAMPLEID}.R1.trimmed.fastq.gz ]
then
echo  fastp step failed with sample ${SAMPLEID}
exit 1
fi

echo `date`
echo stage out trim results of sample ${SAMPLEID}

cp ${SAMPLEID}*trimmed.fastq.gz $OUTPUTDIR/trim
cp ${SAMPLEID}*fastp.* $OUTPUTDIR/QC

echo `date`
echo start ALIGN-SORT with sample ${SAMPLEID}
module load sentieon/202308

echo load bwa index files into memory

#sentieon bwa shm $REFERENCE

echo run alignment step

(sentieon bwa mem -R "@RG\tID:${ID}\tSM:${SM}\tPL:${PL}" \
  -t $THREADS $REFERENCE ${SAMPLEID}.R1.trimmed.fastq.gz ${SAMPLEID}.R2.trimmed.fastq.gz || echo -n 'error' ) \
  | sentieon util sort -r $REFERENCE -o ${SAMPLEID}_SORTED.BAM -t $THREADS --sam2bam -i -
  
  
echo `date`
echo end ALIGN-SORT with sample ${SAMPLEID}

if [ ! -s ${SAMPLEID}_SORTED.BAM ]
then
echo  sentieon bwa failed with sample ${SAMPLEID}
exit 1
fi

echo `date`
echo stage out align results of sample ${SAMPLEID}

cp ${SAMPLEID}_SORTED.* $OUTPUTDIR/sorted

echo `date`
echo start DEDUP with sample ${SAMPLEID}


sentieon driver -t $THREADS -i ${SAMPLEID}_SORTED.BAM \
  --algo LocusCollector --fun score_info ${SAMPLEID}_SCORE.gz

sentieon driver -t $THREADS -i ${SAMPLEID}_SORTED.BAM \
  --algo Dedup --rmdup --score_info ${SAMPLEID}_SCORE.gz --metrics ${SAMPLEID}_DEDUP_METRIC.TXT ${SAMPLEID}_DEDUPED.BAM

echo `date`
echo end DEDUP with sample ${SAMPLEID}

if [ ! -s ${SAMPLEID}_DEDUPED.BAM ]
then
echo Dedup failed with sample ${SAMPLEID}
exit 1
fi

echo `date`
echo stage out dedup results of sample ${SAMPLEID}

cp ${SAMPLEID}_DEDUPED.* $OUTPUTDIR/dedup

echo `date`
echo start VarCalling with sample ${SAMPLEID} to produce output in gVCF format

## skip DNAscope step because we dont have golden variants and indels to build the model with

sentieon driver -t $THREADS -r $REFERENCE \
  -i ${SAMPLEID}_DEDUPED.BAM --algo QualCal  ${SAMPLEID}_RECAL_DATA.TABLE

if [ -z "$DBSNP" ]
then

echo $DBSNP dbSNP is missing. Will call variants without dbSNP parameter

sentieon driver -t $THREADS -r $REFERENCE  -i ${SAMPLEID}_DEDUPED.BAM \
  -q ${SAMPLEID}_RECAL_DATA.TABLE --algo Haplotyper --emit_mode gvcf  \
   ${SAMPLEID}_VARIANT.GVCF

else

echo $DBSNP dbSNP file provided. Will call variants with dbSNP parameter

sentieon driver -t $THREADS -r $REFERENCE  -i ${SAMPLEID}_DEDUPED.BAM \
  -q ${SAMPLEID}_RECAL_DATA.TABLE --algo Haplotyper --emit_mode gvcf   \
  -d $DBSNP  ${SAMPLEID}_VARIANT.GVCF  
fi

  
echo `date`
echo end VarCalling with sample ${SAMPLEID}

if [ ! -s ${SAMPLEID}_VARIANT.GVCF ]
then
echo VarCalling failed with sample ${SAMPLEID}
exit 1
fi

echo `date`
echo stage out variant results of sample ${SAMPLEID}

cp ${SAMPLEID}_VARIANT.* $OUTPUTDIR/Variants

echo `date`
echo start QC of sample ${SAMPLEID}
    
sentieon driver -r $REFERENCE -t $THREADS -i ${SAMPLEID}_SORTED.BAM \
    --algo MeanQualityByCycle ${SAMPLEID}_SORTED_mq_metrics.txt \
    --algo QualDistribution ${SAMPLEID}_SORTED_qd_metrics.txt \
    --algo GCBias --summary ${SAMPLEID}_SORTED_gc_summary.txt ${SAMPLEID}_SORTED_gc_metrics.txt \
    --algo AlignmentStat  --adapter_seq '' ${SAMPLEID}_SORTED_aln_metrics.txt \
    --algo InsertSizeMetricAlgo ${SAMPLEID}_SORTED_is_metrics.txt || { echo "${SAMPLEID}_SORTED Metrics failed"; exit 1; }
    
sentieon driver -r $REFERENCE -t $THREADS -i ${SAMPLEID}_DEDUPED.BAM \
    --algo MeanQualityByCycle ${SAMPLEID}_DEDUPED_mq_metrics.txt \
    --algo QualDistribution ${SAMPLEID}_DEDUPED_qd_metrics.txt \
    --algo GCBias --summary ${SAMPLEID}_DEDUPED_gc_summary.txt ${SAMPLEID}_DEDUPED_gc_metrics.txt \
    --algo AlignmentStat  --adapter_seq '' ${SAMPLEID}_DEDUPED_aln_metrics.txt \
    --algo InsertSizeMetricAlgo ${SAMPLEID}_DEDUPED_is_metrics.txt || { echo "${SAMPLEID}_DEDUPED Metrics failed"; exit 1; }


sentieon plot GCBias -o ${SAMPLEID}_SORTED_gc-report.pdf ${SAMPLEID}_SORTED_gc_metrics.txt
sentieon plot QualDistribution -o ${SAMPLEID}_SORTED_qd-report.pdf ${SAMPLEID}_SORTED_qd_metrics.txt
sentieon plot MeanQualityByCycle -o ${SAMPLEID}_SORTED_mq-report.pdf ${SAMPLEID}_SORTED_mq_metrics.txt
sentieon plot InsertSizeMetricAlgo -o ${SAMPLEID}_SORTED_is-report.pdf ${SAMPLEID}_SORTED_is_metrics.txt


sentieon plot GCBias -o ${SAMPLEID}_DEDUPED_gc-report.pdf ${SAMPLEID}_DEDUPED_gc_metrics.txt
sentieon plot QualDistribution -o ${SAMPLEID}_DEDUPED_qd-report.pdf ${SAMPLEID}_DEDUPED_qd_metrics.txt
sentieon plot MeanQualityByCycle -o ${SAMPLEID}_DEDUPED_mq-report.pdf ${SAMPLEID}_DEDUPED_mq_metrics.txt
sentieon plot InsertSizeMetricAlgo -o ${SAMPLEID}_DEDUPED_is-report.pdf ${SAMPLEID}_DEDUPED_is_metrics.txt

module unload sentieon/202112.04
module load SAMtools/1.17-IGB-gcc-8.2.0
module load BCFtools/1.17-IGB-gcc-8.2.0
module load MultiQC/1.14-IGB-gcc-8.2.0-Python-3.7.2

samtools flagstat ${SAMPLEID}_SORTED.BAM > ${SAMPLEID}_SORTED.BAM.flagstat.stats

samtools idxstats ${SAMPLEID}_SORTED.BAM > ${SAMPLEID}_SORTED.BAM.idxstats.stats

samtools stats ${SAMPLEID}_SORTED.BAM > ${SAMPLEID}_SORTED.BAM.stats

samtools flagstat ${SAMPLEID}_DEDUPED.BAM > ${SAMPLEID}_DEDUPED.BAM.flagstat.stats

samtools idxstats ${SAMPLEID}_DEDUPED.BAM > ${SAMPLEID}_DEDUPED.BAM.idxstats.stats

samtools stats ${SAMPLEID}_DEDUPED.BAM > ${SAMPLEID}_DEDUPED.BAM.stats

bcftools stats -F ${REFERENCE}  ${SAMPLEID}_VARIANT.GVCF > ${SAMPLEID}_VARIANT.GVCF.stats


mkdir ${SAMPLEID}_multiqc

multiqc --outdir ${SAMPLEID}_multiqc . 

echo `date`
echo end QC of sample ${SAMPLEID}

echo `date`
echo stage out QC results of sample ${SAMPLEID}

cp ${SAMPLEID}*txt $OUTPUTDIR/QC

cp ${SAMPLEID}*pdf $OUTPUTDIR/QC

cp ${SAMPLEID}*stats $OUTPUTDIR/QC

cp -r ${SAMPLEID}_multiqc $OUTPUTDIR/QC

echo `date`
echo DONE with analysis of sample ${SAMPLEID}. Reset and exit

unlink $SAMPLER1
unlink $SAMPLER2

cd $OUTPUTDIR
rm -rf $SCRATCH

