#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J cmac_ccs_stats
#SBATCH -C usage_mail

ml load bioinfo-tools
ml load SeqKit/0.15.0

# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}

# Copy the HiFi reads that made the assembly
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/ccsreads/pt_036_001/pt_036_001.ccsreads.fastq.gz ${SNIC_TMP}/

echo ""
ls -lh ${SNIC_TMP}/*
echo ""

seqkit stats ${SNIC_TMP}/pt_036_001.ccsreads.fastq.gz > ${SNIC_TMP}/pt_036_001.ccsreads.stats
cat ${SNIC_TMP}/pt_036_001.ccsreads.stats

# Copy output back to HOME
cp ${SNIC_TMP}/pt_036_001.ccsreads.stats ${HOME}/projects/cmac_XY/DiscoverY/

