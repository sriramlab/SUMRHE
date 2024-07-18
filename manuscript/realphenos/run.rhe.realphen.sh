#!/bin/sh


i=$SGE_TASK_ID

phen=`awk -v idx=$i 'NR==idx' residual_pheno_list.txt`

DIR=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/291k_ukbb_phenotype

gen=/u/project/sgss/UKBB/data/cal/filter4_no_mhc/filter4_no_mhc
src=/u/home/b/bronsonj/project-sriram/RHE-mc/build_sums
out=${DIR}/residual_ldsc/rheouts
pheno=/u/scratch/z/zhengton/SUM_RHE/data/real_data/pheno_res/${phen}.pheno
annot=${DIR}/singleannot.txt

$src/RHEmc_mem -g $gen -p $pheno -k 100 -jn 1000 -o $out/${phen}.rheout.txt -annot $annot
