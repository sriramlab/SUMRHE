#!/bin/bash

folder_path=/u/home/z/zhengton/project-sriram/SUM_RHE/
code_folder=/u/home/z/zhengton/project-sriram/softwares/

## step 1: make grm
k=$1

mkdir -p ${folder_path}/data/grm/
geno=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_50k/geno/50k_split
out_folder=/u/home/z/zhengton/scratch/SUM_RHE/outputs/runtime/greml/
mkdir -p ${out_folder}

${code_folder}/gcta/gcta64 --bfile ${geno} --make-grm --out ${out_folder}/50k_split.${k} --thread-num 6


# step 2: estimate variance
grm=${folder_path}/data/grm/10k_split
for c in sims_0.25_test sims_0.25_test_0.01 sims_0.25_test_0.1 sims_0.4_test sims_0.4_test_0.01 sims_0.4_test_0.1 sims_null_test; do
    sim_folder=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_10k/sims/${c}/
    
    out=/u/home/z/zhengton/scratch/SUM_RHE/outputs/${c}/greml/${k}/
    mkdir -p ${out}
    ${code_folder}/gcta/gcta64 --grm ${grm} --pheno ${sim_folder}/sim_10k_${k}.phen --reml --out ${out}/10k_split --thread-num 10

done