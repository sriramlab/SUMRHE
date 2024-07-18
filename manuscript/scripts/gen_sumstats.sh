#!/bin/sh

source /u/local/Modules/default/init/modules.sh

h2=0.1_test
pol=_0.01
part=_p0.05

suffix=${h2}${pol}${part}
size=10k

i=`expr $SGE_TASK_ID - 1`


psrc=/u/home/b/bronsonj/project-sriram
dir=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${size}

all_bfile=${dir}/geno/${size}_split
phen_dir=${dir}/sims/sims_${suffix}

out_dir=${dir}/sums/sums_${suffix}


echo $out_dir

pheno_file=${phen_dir}/sim_${size}_${i}.phen

mkdir -p ${out_dir}

ldscout=${out_dir}_ldsc

ls ${ldscout}/sim_${size}_${i}.sumstat > /dev/null 2>&1

if [ $? -ne 0 ]
then
    # run GWAS
    ${psrc}/plink2 --bfile ${all_bfile} \
      --pheno ${pheno_file} \
      --variance-standardize \
      --linear allow-no-covars \
      --threads 6 \
      --out ${out_dir}/sim_${size}_${i}.sumstat

    sed -i '1s/^#//' ${out_dir}/sim_${size}_${i}.sumstat.pheno.glm.linear

    # convert to ldsc format    
    module load R
    mkdir -p ${ldscout}

    Rscript ${dir}/format_ldsc_sumstat2.R ${out_dir}/sim_${size}_${i}.sumstat.pheno.glm.linear ${ldscout}/sim_${size}_${i}.sumstat
    rm ${out_dir}/sim_${size}_${i}.sumstat.pheno.glm.linear
fi

