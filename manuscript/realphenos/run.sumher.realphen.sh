#!/bin/bash


folder_path=/u/home/z/zhengton/project-sriram/SUM_RHE/
ldak=/u/home/z/zhengton/project-sriram/softwares/ldak5.2.linux
plink2=/u/home/z/zhengton/project-sriram/softwares/plink2

geno=/u/home/z/zhengton/project-ukbiobank/data/geno/cal/filter4_no_mhc/filter4_no_mhc
pheno_name=$1
pheno=/u/project/sriram/alipazok/Data/new_ukkb_phenotypes/${pheno_name}.pheno
covar=/u/project/sriram/alipazok/Data/new_ukkb_phenotypes/${pheno_name}.covar
k=$1

ref=/u/home/z/zhengton/scratch/SUM_RHE/outputs/real_data/sumher/${pheno_name}/
out=/u/home/z/zhengton/scratch/SUM_RHE/outputs/real_data/sumher_gcta/${pheno_name}/
mkdir -p ${out}

${plink2} --pheno-quantile-normalize --covar-variance-standardize --pheno ${pheno} --covar ${covar} --bfile ${geno} --glm 'hide-covar' --out ${out}/quant --threads 6

f=${out}/quant.pheno.glm.linear
sed -i '1s/REF/A2/' "$f"
sed -i '1s/OBS_CT/n/' "$f"
sed -i '1s/T_STAT/Z/' "$f"
sed -i '1s/ID/Predictor/' "$f"

${ldak} --sum-hers ${out}/snpher --summary ${ref}/quant.pheno.glm.linear --tagfile ${folder_path}/data/sumher/10k_split_GCTA-Thin.tagging