#!/bin/bash -l
 
#SBATCH -A snic2021-5-153
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J my_job
#SBATCH -C usage_mail

# TMP Directory
#$SNIC_TMP

input=$1

echo "${input}_1 ${input}_2"
