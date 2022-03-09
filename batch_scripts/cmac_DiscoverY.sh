#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 3-00:00:00
#SBATCH -J cmac_DiscoverY
#SBATCH -C mem1TB
#SBATCH -C usage_mail


# TMP Directory
#$SNIC_TMP

#JOB ID
#$SLURM_JOB_ID

# Move to the SNIC_TMP directory
cd ${SNIC_TMP}


# Copy the DiscoverY tools over
cp -r ${HOME}/packages/DiscoverY ${SNIC_TMP}/
rm -r ${SNIC_TMP}/DiscoverY/data
mkdir ${SNIC_TMP}/DiscoverY/data
echo ""
ls -lh ${SNIC_TMP}/*
ls -lh ${SNIC_TMP}/DiscoverY/data/
echo ""


# Copy the genome assembly
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/DiscoverY/data/male_contigs.fasta
echo ""
ls -lh ${SNIC_TMP}/*
ls -lh ${SNIC_TMP}/DiscoverY/data/
echo ""


# Copy the HiFi reads that made the assembly
#cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/ccsreads/pt_036_001/pt_036_001.ccsreads.fastq.gz ${SNIC_TMP}/male_reads.fastq.gz
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/rawdata/pt_036_001/r64077_20201013_114036/2_B01/m64077_201014_183827.subreads.fq.gz ${SNIC_TMP}/male_reads.fastq.gz
gunzip ${SNIC_TMP}/male_reads.fastq.gz
echo ""
ls -lh ${SNIC_TMP}/*
ls -lh ${SNIC_TMP}/DiscoverY/data/
echo ""


# Make K-mers from HiFi reads
cd ${SNIC_TMP}/DiscoverY/dependency

./run_dsk_Linux.sh ${SNIC_TMP}/male_reads.fastq 25
echo ""
ls -lh ${SNIC_TMP}/*
ls -lh ${SNIC_TMP}/DiscoverY/data/
echo ""

cp ${SNIC_TMP}/DiscoverY/dependency/dsk_output/kmers_from_reads ${SNIC_TMP}/DiscoverY/data/kmers_from_male_reads


# Copy the files for female F1 sequenced in Sayadi et al. 2019 and concatenate
cp /proj/snic2021-6-30/sayadi_et_al_2019_illumina_reads/CmacF1*.fastq.gz ${SNIC_TMP}/

zcat ${SNIC_TMP}/CmacF1*.fastq.gz > ${SNIC_TMP}/female_reads.fastq
echo ""
ls -lh ${SNIC_TMP}/*
ls -lh ${SNIC_TMP}/DiscoverY/data/
echo ""

rm ${SNIC_TMP}/CmacF1*.fastq.gz
echo ""
ls -lh ${SNIC_TMP}/*
ls -lh ${SNIC_TMP}/DiscoverY/data/
echo ""



# Make K-mers from Female Illumina reads
cd ${SNIC_TMP}/DiscoverY/dependency
rm -r ${SNIC_TMP}/DiscoverY/dependency/dsk_output

./run_dsk_Linux.sh ${SNIC_TMP}/female_reads.fastq 25
echo ""
ls -lh ${SNIC_TMP}/*
ls -lh ${SNIC_TMP}/DiscoverY/data/
echo ""

cp ${SNIC_TMP}/DiscoverY/dependency/dsk_output/kmers_from_reads ${SNIC_TMP}/DiscoverY/data/female_kmers



# Run DiscoverY
cd ${SNIC_TMP}/DiscoverY
echo "Running DiscoverY..."
python3 discoverY.py --kmer_size 25 --mode female+male --female_kmers_set --female_bloom_capacity 1000000000
echo ""
ls -lh ${SNIC_TMP}/*
ls -lh ${SNIC_TMP}/DiscoverY/data/
echo ""



# Copy output back to HOME
cp ${SNIC_TMP}/DiscoverY/proportion_annotated_contigs.fasta ${HOME}/projects/cmac_XY/DiscoverY/
#cp ${SNIC_TMP}/DiscoverY/data/* ${HOME}/projects/cmac_XY/DiscoverY/data/
