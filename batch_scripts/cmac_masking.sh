#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH -J cmac_masking
#SBATCH -C usage_mail


# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

ml load bwa/0.7.17
ml load samtools/1.14
ml load RepeatMasker/4.1.2-p1

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}

# Copy the genome assembly
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/

# Copy the repeat library
cp /home/axelw/data/galbraith_repeats/rm1.2_beetle_library.lib ${SNIC_TMP}/

## Masking repeats
echo "Masking repeats..."
mkdir ${SNIC_TMP}/rm_results

RepeatMasker -pa 5 -e rmblast -xsmall -gff -dir ${SNIC_TMP}/rm_results -lib ${SNIC_TMP}/rm1.2_beetle_library.lib ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta

echo "Masking done."

ls -lh ${SNIC_TMP}/
ls -lh ${SNIC_TMP}/rm_results

# Copy results back
cp -r ${SNIC_TMP}/rm_results /home/axelw/data/galbraith_repeats/

# Cleanup
rm -r ${SNIC_TMP}/*
echo "Done"








