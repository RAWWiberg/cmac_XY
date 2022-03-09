#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 5-00:00:00
#SBATCH -J cmac_F_filtering
#SBATCH -C usage_mail


# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

ml load bwa/0.7.17
ml load samtools/1.14

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}

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
# skip alignmenst with MAPQ < 1: -q 1
# keep properly paired: -f 2
# exclude alignments that are not primary: -F 2048
samtools view -b -@ 20 -f 2 -F 2048 -q 1 ${SNIC_TMP}/CmacF1.bam > ${SNIC_TMP}/CmacF1_filt.bam
samtools index -@ 20 ${SNIC_TMP}/CmacF1_filt.bam
samtools idxstats ${SNIC_TMP}/CmacF1_filt.bam > ${SNIC_TMP}/CmacF1_filt.idxstats
samtools flagstat ${SNIC_TMP}/CmacF1_filt.bam > ${SNIC_TMP}/CmacF1_filt.flagstats

## Copy output back to home
cp ${SNIC_TMP}/CmacF1_filt.idxstats ${HOME}/projects/cmac_XY/SATC/
cp ${SNIC_TMP}/CmacF1_filt.flagstats ${HOME}/projects/cmac_XY/SATC/

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
# skip alignmenst with MAPQ < 1: -q 1
# keep properly paired: -f 2
# exclude supplementary alignments: -F 2048 
samtools view -b -@ 20 -f 2 -F 2048 -q 1 ${SNIC_TMP}/CmacF2.bam > ${SNIC_TMP}/CmacF2_filt.bam
samtools index -@ 20 ${SNIC_TMP}/CmacF2_filt.bam
samtools idxstats ${SNIC_TMP}/CmacF2_filt.bam > ${SNIC_TMP}/CmacF2_filt.idxstats
samtools flagstat ${SNIC_TMP}/CmacF2_filt.bam > ${SNIC_TMP}/CmacF2_filt.flagstats

## Copy output back to home
cp ${SNIC_TMP}/CmacF2_filt.idxstats ${HOME}/projects/cmac_XY/SATC/
cp ${SNIC_TMP}/CmacF2_filt.flagstats ${HOME}/projects/cmac_XY/SATC/

## Cleanup
rm ${SNIC_TMP}/CmacF2*
echo "Done"








