
source /u/local/Modules/default/init/modules.sh
module load anaconda2/2019.10


path=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/291k_ukbb_phenotype

while IFS= read -r line; do
    phen=$line
    echo $phen

    DIR=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/291k_ukbb_phenotype
    ld_dir=${DIR}/ldscores_281k
    src_dir=/u/home/b/bronsonj/project-sriram/ldsc
    out=${path}/ldscouts_twostep/${phen}_ldscout

    # run ldsc
    sumstat=${path}/sumstats_ldsc/${phen}_ldsc.sumstat

    python "$src_dir"/ldsc.py \
          --h2 $sumstat\
          --ref-ld ${ld_dir}/all \
          --w-ld ${ld_dir}/all \
          --out $out \
          --not-M-5-50 \
          --n-blocks 1000 \
          --chisq-max 99999 \
          --two-step 99999

done < "$path/residual_pheno_list.txt"
