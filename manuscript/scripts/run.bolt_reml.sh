#!/bin/sh
ldscores=/path/to/project/SUM_RHE/data/all.l2.ldscore.gz 
src=/path/to/project/BOLT-LMM_v2.3.5
size=50k
dir=/path/to/project/summaryRHE/300k/split_${size}
i=$1
phenoCol=pheno
folder_path=/path/to/project/SUM_RHE/

bfile=${dir}/geno/${size}_split
model_snp_file=${folder_path}/data/bolt/modelSnps.txt

for c in "sims_0.25_test" "sims_0.4_test" "sims_0.25_test_0.1_p0.05" "sims_0.4_test_0.1_p0.05" "sims_0.25_test_0.01_p0.05" "sims_0.4_test_0.01_p0.05" "sims_null_test"; do
    outpath=/path/to/scratch/SUM_RHE/outputs/${c}/bolt-reml
    mkdir -p ${outpath}/${i}
    phen=${dir}/sims/${c}

    $src/bolt --bfile=$bfile --phenoFile=$phen/sim_10k_${i}.phen --phenoCol=$phenoCol --lmm --LDscoresFile=$ldscores --statsFile=${outpath}/${i}/stats.tab \
        2>&1 | tee $outpath/${i}/out.log
    time $src/bolt --bfile=$bfile --phenoFile=$phen/sim_50k_${i}.phen --phenoCol=$phenoCol --reml --modelSnps=$model_snp_file --numThreads=6 \
        2>&1 | tee $outpath/${i}/out.log
done