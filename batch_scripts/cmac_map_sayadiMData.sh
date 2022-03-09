#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH -J cmac_Mmapping
#SBATCH -C usage_mail


# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

ml load bwa/0.7.17
ml load samtools/1.14

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}

# Copy the genome assembly
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/

# Make bwa index
bwa index ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta


# M1
echo "Mapping for M1..."
# Copy the Illumina files and concatenate.
echo "Copying"
echo ""
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacM1*.fastq.gz ${SNIC_TMP}/
ls -lh ${SNIC_TMP}/

echo "Concatenating"
echo ""
zcat ${SNIC_TMP}/CmacM1*_R1_*.fastq.gz > ${SNIC_TMP}/CmacM1_R1.fq
zcat ${SNIC_TMP}/CmacM1*_R2_*.fastq.gz > ${SNIC_TMP}/CmacM1_R2.fq
rm *fastq.gz
ls -lh ${SNIC_TMP}/

# Map reads
echo "Mapping"
echo ""
bwa mem -t 10 ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/CmacM1_R1.fq ${SNIC_TMP}/CmacM1_R2.fq | samtools sort - | samtools view - -b > ${SNIC_TMP}/CmacM1.bam
samtools index ${SNIC_TMP}/CmacM1.bam
samtools idxstats ${SNIC_TMP}/CmacM1.bam > ${SNIC_TMP}/CmacM1.idxstats
samtools flagstat ${SNIC_TMP}/CmacM1.bam > ${SNIC_TMP}/CmacM1.flagstats

# Copy bam files back to /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads
cp ${SNIC_TMP}/CmacM1.bam /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/
cp ${SNIC_TMP}/CmacM1.bam.bai /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/

# Copy output back to home
cp ${SNIC_TMP}/CmacM1.idxstats ${HOME}/projects/cmac_XY/SATC/
cp ${SNIC_TMP}/CmacM1.flagstats ${HOME}/projects/cmac_XY/SATC/

## Cleanup
rm ${SNIC_TMP}/CmacM1*
echo "Done"



#F2
echo "Mapping for M2..."
# Copy the Illumina files and concatenate.
echo "Copying"
echo ""
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacM2*.fastq.gz ${SNIC_TMP}/
echo "Concatenating"
echo ""
zcat ${SNIC_TMP}/CmacM2*_R1_*.fastq.gz > ${SNIC_TMP}/CmacM2_R1.fq
zcat ${SNIC_TMP}/CmacM2*_R2_*.fastq.gz > ${SNIC_TMP}/CmacM2_R2.fq
rm *fastq.gz
ls -lh ${SNIC_TMP}/

# Map reads
echo "Mapping"
echo ""
bwa mem -t 20 ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/CmacM2_R1.fq ${SNIC_TMP}/CmacM2_R2.fq | samtools sort -@ 20 - | samtools view -@ 20 - -b > ${SNIC_TMP}/CmacM2.bam
samtools index ${SNIC_TMP}/CmacM2.bam
samtools idxstats ${SNIC_TMP}/CmacM2.bam > ${SNIC_TMP}/CmacM2.idxstats
samtools flagstat ${SNIC_TMP}/CmacM2.bam > ${SNIC_TMP}/CmacM2.flagstats

# Copy bam files back to /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads
cp ${SNIC_TMP}/CmacM2.bam /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/
cp ${SNIC_TMP}/CmacM2.bam.bai /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/

# Copy output back to home
cp ${SNIC_TMP}/CmacM2.idxstats ${HOME}/projects/cmac_XY/SATC/
cp ${SNIC_TMP}/CmacM2.flagstats ${HOME}/projects/cmac_XY/SATC/

# Cleanup
rm ${SNIC_TMP}/CmacM2*
echo "Done"








