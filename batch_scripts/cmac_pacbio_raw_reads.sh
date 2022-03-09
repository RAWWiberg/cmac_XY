#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J cmac_raw_reads
#SBATCH -C usage_mail

ml load bioinfo-tools
ml load SeqKit/0.15.0
ml load samtools/1.14

# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}

# Copy the HiFi reads that made the assembly
#cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/rawdata/pt_036_001/r64077_20201013_114036/2_B01/m64077_201014_183827.subreads.bam ${SNIC_TMP}/
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/rawdata/pt_036_001/r64077_20201013_114036/2_B01/m64077_201014_183827.subreads.fq.gz ${SNIC_TMP}/

echo ""
ls -lh ${SNIC_TMP}/*
echo ""

# Convert .bam file to fastq
#samtools fastq -0 ${SNIC_TMP}/m64077_201014_183827.subreads.fq ${SNIC_TMP}/m64077_201014_183827.subreads.bam
#gzip ${SNIC_TMP}/m64077_201014_183827.subreads.fq

echo ""
ls -lh ${SNIC_TMP}/*
echo ""

# Collect stats
seqkit stats ${SNIC_TMP}/m64077_201014_183827.subreads.fq.gz > ${SNIC_TMP}/m64077_201014_183827.subreads.stats
cat ${SNIC_TMP}/m64077_201014_183827.subreads.stats

# Copy output back to HOME
cp ${SNIC_TMP}/m64077_201014_183827.subreads.stats ${HOME}/projects/cmac_XY/DiscoverY/

#cp ${SNIC_TMP}/m64077_201014_183827.subreads.fq.gz /proj/snic2021-6-30/delivery04381/INBOX/pt_036/rawdata/pt_036_001/r64077_20201013_114036/2_B01/



