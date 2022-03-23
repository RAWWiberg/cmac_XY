#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 5-00:00:00
#SBATCH -J CmacM_filtering_pt_036_002
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
cp /home/axelw/data/galbraith_repeats/rm_results_pt_036_002/pt_036_002.hifiasm20210112.primary.fasta.out.bed ${SNIC_TMP}/


## CmacM1
echo "Filtering CmacM1..."
# Copy the bam files.
echo "Copying"
echo ""
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacM1_pt_036_002.bam ${SNIC_TMP}/
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacM1_pt_036_002.bam.bai ${SNIC_TMP}/
ls -lh ${SNIC_TMP}/

# Filter
echo "Filtering"
echo ""
# remove reads that occur in repeatmasked regions
bedtools intersect -abam ${SNIC_TMP}/CmacM1_pt_036_002.bam -b ${SNIC_TMP}/pt_036_002.hifiasm20210112.primary.fasta.out.bed -v > ${SNIC_TMP}/CmacM1_pt_036_002_filtreps.bam

samtools view -b -@ 10 -f 2 -F 2048 -q 1 ${SNIC_TMP}/CmacM1_pt_036_002_filtreps.bam > ${SNIC_TMP}/CmacM1_pt_036_002_filtrepsreads.bam

samtools index -@ 10 ${SNIC_TMP}/CmacM1_pt_036_002_filtrepsreads.bam

samtools idxstats ${SNIC_TMP}/CmacM1_pt_036_002_filtrepsreads.bam > ${SNIC_TMP}/CmacM1_pt_036_002_filtrepsreads.idxstats

samtools flagstat ${SNIC_TMP}/CmacM1_pt_036_002_filtrepsreads.bam > ${SNIC_TMP}/CmacM1_pt_036_002_filtrepsreads.flagstats

## Copy output back to home
cp ${SNIC_TMP}/CmacM1_pt_036_002_filtrepsreads.idxstats ${HOME}/projects/XY/SATC/
cp ${SNIC_TMP}/CmacM1_pt_036_002_filtrepsreads.flagstats ${HOME}/projects/XY/SATC/

## Cleanup
rm ${SNIC_TMP}/CmacM1*
echo "Done"



## CmacM2
echo "Filtering CmacM2..."
# Copy the bam files.
echo "Copying"
echo ""
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacM2_pt_036_002.bam ${SNIC_TMP}/
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacM2_pt_036_002.bam.bai ${SNIC_TMP}/
ls -lh ${SNIC_TMP}/

# Filter
echo "Filtering"
echo ""
# remove reads that occur in repeatmasked regions
bedtools intersect -abam ${SNIC_TMP}/CmacM2_pt_036_002.bam -b ${SNIC_TMP}/pt_036_002.hifiasm20210112.primary.fasta.out.bed -v > ${SNIC_TMP}/CmacM2_pt_036_002_filtreps.bam

samtools view -b -@ 10 -f 2 -F 2048 -q 1 ${SNIC_TMP}/CmacM2_pt_036_002_filtreps.bam > ${SNIC_TMP}/CmacM2_pt_036_002_filtrepsreads.bam

samtools index -@ 10 ${SNIC_TMP}/CmacM2_pt_036_002_filtrepsreads.bam

samtools idxstats ${SNIC_TMP}/CmacM2_pt_036_002_filtrepsreads.bam > ${SNIC_TMP}/CmacM2_pt_036_002_filtrepsreads.idxstats

samtools flagstat ${SNIC_TMP}/CmacM2_pt_036_002_filtrepsreads.bam > ${SNIC_TMP}/CmacM2_pt_036_002_filtrepsreads.flagstats

## Copy output back to home
cp ${SNIC_TMP}/CmacM2_pt_036_002_filtrepsreads.idxstats ${HOME}/projects/XY/SATC/
cp ${SNIC_TMP}/CmacM2_pt_036_002_filtrepsreads.flagstats ${HOME}/projects/XY/SATC/

## Cleanup
rm ${SNIC_TMP}/CmacM2*
echo "Done"








