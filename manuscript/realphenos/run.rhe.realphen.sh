#!/bin/sh


i=$SGE_TASK_ID

phen=`awk -v idx=$i 'NR==idx' residual_pheno_list.txt`

DIR=/path/to/project/summaryRHE/300k/291k_ukbb_phenotype

gen=/path/to/authorized_data/data/cal/filter4_no_mhc/filter4_no_mhc
src=/path/to/project/RHE-mc/build_sums
out=${DIR}/residual_ldsc/rheouts
pheno=/path/to/scratch/SUM_RHE/data/real_data/pheno_res/${phen}.pheno
annot=${DIR}/singleannot.txt

$src/RHEmc_mem -g $gen -p $pheno -k 100 -jn 1000 -o $out/${phen}.rheout.txt -annot $annot
