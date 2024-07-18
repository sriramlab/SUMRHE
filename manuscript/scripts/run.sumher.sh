#!/bin/bash

k=$1
for sample_size in "10k" "50k"; do
geno=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${sample_size}/geno/${sample_size}_split
    for c in sims_0.25_test_0.01 sims_0.25_test_0.1 sims_0.25_test sims_0.4_test_0.01 sims_0.4_test_0.1 sims_0.4_test sims_null_test; do
        sim_folder=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${sample_size}/sims/${c}/

        folder_path=/u/home/z/zhengton/project-sriram/SUM_RHE/
        ldak=/u/home/z/zhengton/project-sriram/softwares/ldak5.2.linux


        out=/u/home/z/zhengton/scratch/SUM_RHE/outputs/${c}/sumher_gcta/${k}/
        if [[ $sample_size == "50k" ]]; then
            out=/u/home/z/zhengton/scratch/SUM_RHE/outputs/50k/${c}/sumher_gcta/${k}/
        fi

        mkdir -p ${out}

        ${ldak} --linear ${out}/quant --bfile ${geno} --pheno ${sim_folder}/sim_${sample_size}_${k}.phen --max-threads 10

        ${ldak} --sum-hers ${out}/snpher --summary ${out}/quant.summaries --tagfile ${folder_path}/data/sumher/${sample_size}_split_GCTA-Thin.tagging --max-threads 10
        rm ${out}/quant.*
    done
done