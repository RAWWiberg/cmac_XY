#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 5-00:00:00
#SBATCH -J cmac_F_filtering
#SBATCH -C usage_mail


# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

ml load samtools/1.14
ml load BEDTools/2.29.2

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}

# Copy over the repeats file
cp /home/axelw/data/galbraith_repeats/rm_results/pt_036_001_hifiasm_20201223.primary.fasta.out.bed ${SNIC_TMP}/


## F1
echo "Filtering F1..."
# Copy the bam files.
echo "Copying"
echo ""
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacF1.bam ${SNIC_TMP}/
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacF1.bam.bai ${SNIC_TMP}/
ls -lh ${SNIC_TMP}/

# Filter
echo "Filtering"
echo ""
# remove reads that occur in repeatmasked regions
bedtools intersect -abam ${SNIC_TMP}/CmacF1.bam -b ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta.out.bed -v > ${SNIC_TMP}/CmacF1_filtreps.bam
samtools view -b -@ 10 -f 2 -F 2048 -q 1 ${SNIC_TMP}/CmacF1_filtreps.bam > ${SNIC_TMP}/CmacF1_filtrepsreads.bam
samtools index -@ 10 ${SNIC_TMP}/CmacF1_filtrepsreads.bam
samtools idxstats ${SNIC_TMP}/CmacF1_filtrepsreads.bam > ${SNIC_TMP}/CmacF1_filtrepsreads.idxstats
samtools flagstat ${SNIC_TMP}/CmacF1_filtrepsreads.bam > ${SNIC_TMP}/CmacF1_filtrepsreads.flagstats

## Copy output back to home
cp ${SNIC_TMP}/CmacF1_filtrepsreads.idxstats ${HOME}/projects/cmac_XY/SATC/
cp ${SNIC_TMP}/CmacF1_filtrepsreads.flagstats ${HOME}/projects/cmac_XY/SATC/

## Cleanup
rm ${SNIC_TMP}/CmacF1*
echo "Done"



## F2
echo "Filtering F2..."
# Copy the bam files.
echo "Copying"
echo ""
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacF2.bam ${SNIC_TMP}/
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacF2.bam.bai ${SNIC_TMP}/
ls -lh ${SNIC_TMP}/

# Filter
echo "Filtering"
echo ""
# remove reads that occur in repeatmasked regions
bedtools intersect -abam ${SNIC_TMP}/CmacF2.bam -b ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta.out.bed -v > ${SNIC_TMP}/CmacF2_filtreps.bam
samtools view -b -@ 10 -f 2 -F 2048 -q 1 ${SNIC_TMP}/CmacF2_filtreps.bam > ${SNIC_TMP}/CmacF2_filtrepsreads.bam
samtools index -@ 10 ${SNIC_TMP}/CmacF2_filtrepsreads.bam
samtools idxstats ${SNIC_TMP}/CmacF2_filtrepsreads.bam > ${SNIC_TMP}/CmacF2_filtrepsreads.idxstats
samtools flagstat ${SNIC_TMP}/CmacF2_filtrepsreads.bam > ${SNIC_TMP}/CmacF2_filtrepsreads.flagstats

## Copy output back to home
cp ${SNIC_TMP}/CmacF2_filtrepsreads.idxstats ${HOME}/projects/cmac_XY/SATC/
cp ${SNIC_TMP}/CmacF2_filtrepsreads.flagstats ${HOME}/projects/cmac_XY/SATC/

## Cleanup
rm ${SNIC_TMP}/CmacF2*
echo "Done"








