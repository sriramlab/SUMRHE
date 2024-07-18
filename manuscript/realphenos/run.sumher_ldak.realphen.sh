#!/bin/bash


folder_path=/u/home/z/zhengton/project-sriram/SUM_RHE/
ldak=/u/home/z/zhengton/project-sriram/softwares/ldak5.2.linux
plink2=/u/home/z/zhengton/project-sriram/softwares/plink2

geno=/u/home/z/zhengton/project-ukbiobank/data/geno/cal/filter4_no_mhc/filter4_no_mhc
pheno_name=$1
pheno=/u/project/sriram/alipazok/Data/new_ukkb_phenotypes/${pheno_name}.pheno
covar=/u/project/sriram/alipazok/Data/new_ukkb_phenotypes/${pheno_name}.covar

out=/u/home/z/zhengton/scratch/SUM_RHE/outputs/real_data/sumher_ldak/${pheno_name}/
mkdir -p ${out}

${plink2} --pheno-quantile-normalize --covar-variance-standardize --pheno ${pheno} --covar ${covar} --bfile ${geno} --glm 'hide-covar' --out ${out}/quant --threads 6

f=${out}/quant.pheno.glm.linear
sed -i '1s/REF/A2/' "$f"
sed -i '1s/OBS_CT/n/' "$f"
sed -i '1s/T_STAT/Z/' "$f"
sed -i '1s/ID/Predictor/' "$f"
${ldak} --thin ${folder_path}/data/sumher/10k_split_thin --bfile ${geno} --window-prune .98 --window-kb 100
awk < ${folder_path}/data/sumher/10k_split_thin.in '{print $1, 1}' > ${folder_path}/data/sumher/10k_split_thin_weights.thin
${ldak} --calc-tagging ${folder_path}/data/sumher/10k_split_LDAK-Thin --bfile ${geno} --weights ${folder_path}/data/sumher/10k_split_thin_weights.thin --power -.25 --window-kb 1000

${ldak} --sum-hers ${out}/snpher --summary ${out}/quant.pheno.glm.linear --tagfile ${folder_path}/data/sumher/10k_split_LDAK-Thin.tagging