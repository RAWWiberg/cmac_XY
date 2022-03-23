### Table of contents
0. [Packages and software](#0packages)
1. [Genome Hi-C scaffolding](#genome_scaffolding)
	1. [Primary assembly](#prime_assembly)
	2. [Purge haplotigs](#purge_hapl)
	3. [Arrow polishing](#arrow)
	4. [Pilon polishing](#pilon)
	5. [Final assembly](#final_assembly)
2. [Sex-chromosome identification](#sex_chrom_id)
	1. [Repeatmasking](#repeatmasking)
	2. [Re-mapping reads to final genome assembly](#re_mapping)
X. [References](#references)



### 0. Packages and software <a name="0packages"></a>
Here are some details on software used


### 1. Genome Hi-C scaffolding <a name="genome_scaffolding"></a>

Hi-C data are from the GA lab. Generated for an inbred line from a different population

Here I am using this data to scaffold the HiFi assemblies generated from the Y-chromosom introgression lines.

#### 1.1 Read mapping

First map the reads. I follow the procedure outlined from [Arima](https://github.com/ArimaGenomics/mapping_pipeline/blob/master/Arima_Mapping_UserGuide_A160156_v02.pdf)

	REF=pt_036_001_hifiasm_20201223.primary.fasta
	PRFX=pt_036_001_hifiasm_20201223.primary
	readsDir=/proj/snic2020-6-128/private/backup/c.mac/raw_data/genome/HiC/P13602/P13602_1001/02-FASTQ/190702_A00689_0045_BHKFFMDSXX/
	readPrfx=P13602_1001_S1_L004

Index the genome

	bwa index -a bwtsw -p ${PRFX} ${REF}
	samtools index ${REF}
	
Align the reads, each read independently, they are paired later

	bwa mem -t 20 ${REF} ${R1}.fq.gz | samtools view -@ 20 -Sb - > ${R1}.bam

	bwa mem -t 20 ${REF} ${R2}.fq.gz | samtools view -@ 20 -Sb - > ${R2}.bam

#### 1.2 Alignment filtering

Some single end mapped reads map to ligation junctions. In these cases, try to identify the regions only mapping to "either side" of the junction. Uses the perl script from Arima `filter_five_end.pl`

	samtools view -h ${R1}.bam | perl filer_five_end.pl | samtools view -Sb > ${R1}_5filt_.bam

	samtools view -h ${R2}.bam | perl filer_five_end.pl | samtools view -Sb > ${R2}_5filt_.bam



#### 1.3 Pair the reads, filter by mapping quality, add readgroup

Reads are paired with the Arima script `two_read_bam_combiner.pl`

	perl two_read_bam_combiner.pl ${R1}_5filt.bam ${R2}_5filt.bam $(which samtools) 10 | samtools view -Sb -t ${REF}.fai - | samtools sort -@ 20 -o ${R}.bam
	
Add read groups with picard

	java -Xmx4G -jar picard.jar AddOrReplaceReadGroups INPUT=${R}.bam OUTPUT=${R}_rdgrp.bam ID={R} LB=${R} SM=${R} PL=ILLUMINA PU=none
	
#### 1.4 Mark PCR duplicates

	java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=./temp/ -jar picard.jar MarkDuplicates INPUT=${R}_rdgrp.bam OUTPUT=${R}_rdgrp_markdup.bam METRICS_FILE=${}

