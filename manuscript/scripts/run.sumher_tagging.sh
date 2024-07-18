#!/bin/bash

folder_path=/u/home/z/zhengton/project-sriram/SUM_RHE/
ldak=/u/home/z/zhengton/project-sriram/softwares/ldak5.2.linux


k=$1
if [[ $k == 1 ]]; then
sample_size="50k"
geno=/u/scratch/b/bronsonj/ref_geno/241k_split
elif [[ $k == 2 ]]; then
sample_size="10k"
geno=/u/scratch/b/bronsonj/ref_geno/281k_split
fi

echo $geno
echo $sample_size

${ldak} --calc-tagging ${folder_path}/data/sumher/${sample_size}_split_GCTA-Thin --bfile ${geno} --power -1 --window-kb 1000 --max-threads 6