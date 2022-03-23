#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH -J cmac_HiCmapping_pt_036_001
#SBATCH -C usage_mail
#SBATCH --mail-user=axel.wiberg@ebc.uu.se
#SBATCH --mail-type=ALL

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
bwa index -a bwtsw -p ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta


## R2
echo "Mapping for R2..."
# Copy the HiC R2 reads.
echo "Copying file over..."
echo ""
cp /proj/snic2020-6-128/private/backup/c.mac/raw_data/genome/HiC/P13602/P13602_1001/02-FASTQ/190702_A00689_0045_BHKFFMDSXX/P13602_1001_S1_L004_R2_001.fastq.gz ${SNIC_TMP}/CmacHiC_R2.fq.gz

ls -lh ${SNIC_TMP}/

# Map reads
echo "Mapping"
echo ""
bwa mem -t 20 ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/CmacHiC_R2.fq.gz | samtools view -@ 20 -Sb - > ${SNIC_TMP}/CmacHiC_R2_pt_036_001.bam
bwa index ${SNIC_TMP}/CmacHiC_R2_pt_036_001.bam
samtools indxstats ${SNIC_TMP}/CmacHiC_R2_pt_036_001.bam > ${SNIC_TMP}/CmacHiC_R2_pt_036_001.indxstats
samtools flagstats ${SNIC_TMP}/CmacHiC_R2_pt_036_001.bam > ${SNIC_TMP}/CmacHiC_R2_pt_036_001.flagstats

# Copy bam files back to /proj/snic2021-6-30/HiC
cp ${SNIC_TMP}/CmacHiC_R2_pt_036_001.bam /proj/snic2021-6-30/HiC_scaffolding

## Copy output back to home
cp ${SNIC_TMP}/CmacHiC_R2_pt_036_001.idxstats ${HOME}/projects/XY/HiC_scaffolding
cp ${SNIC_TMP}/CmacHiC_R2_pt_036_001.flagstats ${HOME}/projects/XY/HiC_scaffolding

echo "Done"








