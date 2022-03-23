#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH -J cmac_Fmapping
#SBATCH -C usage_mail
#SBATCH -C usage_mail
#SBATCH --mail-user=axel.wiberg@ebc.uu.se
#SBATCH --mail-type=ALL

# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

ml load bwa/0.7.17
ml load samtools/1.14


REF=pt_036_001_hifiasm_20201223.primary.fasta
PRFX=pt_036_001_hifiasm_20201223.primary
	readsDir=/proj/snic2020-6-128/private/backup/c.mac/raw_data/genome/HiC/P13602/P13602_1001/02-FASTQ/190702_A00689_0045_BHKFFMDSXX/

readPrfx=P13602_1001_S1_L004

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}

# Copy the genome assembly
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/${REF} ${SNIC_TMP}/

# R1 reads
echo "Mapping R1 reads"

# Copy the HiC reads
cp ${readsDir}/${readPrfx}_R1_001.fastq.gz ${SNIC_TMP}/

# Make bwa index
bwa index -a bwtsw -p ${PRFX} ${SNIC_TMP}/${REF}
samtools index ${SNIC_TMP}/${REF}

# Map the reads

bwa mem -t 20 ${REF} ${readPrfx}_R1_001.fastq.gz | samtools view -@ 20 -Sb - > ${readPrfx}_R1.bam