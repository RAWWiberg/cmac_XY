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

First map the reads. I follow the procedure outlined at []()

