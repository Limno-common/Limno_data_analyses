#!/bin/bash -l
#SBATCH -A snic2021-22-319
#SBATCH -J raxml
#SBATCH -o /proj/snic2020-15-96/private/SlurmOut/raxml.out
#SBATCH -t 5:00:00
#SBATCH -p core
#SBATCH -n 6

module load bioinfo-tools
module load raxml/8.2.12-gcc

raxmlHPC -s /proj/snic2020-15-96/private/RAxML/ChytridSeqsAlign.phy -n ChytridRAxMLtree -m GTRGAMMA -o AY601707.1_Rozella_sp._JEL347_isolate_AFTOL-ID_16,AY635838.1_Rozella_allomycis_isolate_AFTOL-ID_297 -T 6 -f a -N autoMRE -p 12345 -x 12345
