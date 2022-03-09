#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J cmac_masking
#SBATCH -C usage_mail


# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

ml load bioinfo-tools
ml load biopython/1.76-py3

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}

# Copy the genome assembly
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/

## Compute stats
echo "Computing stats..."
python3 ~/bin/fastaStats.py -t p -s n -f ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta > ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta.stats
echo "Done."

ls -lh ${SNIC_TMP}/
ls -lh ${SNIC_TMP}/rm_results

# Copy results back
cp -r ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta.stats /home/axelw/data/

# Cleanup
rm -r ${SNIC_TMP}/*
echo "Done"








