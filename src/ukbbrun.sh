




path=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/291k_ukbb_phenotype/residual_ldsc
dir=/u/home/b/bronsonj/project-sriram/sumrhe
src=${dir}/src
demo=${dir}/data
outdir=${path}/sumrheouts

mkdir -p ${outdir}

while IFS= read -r line; do
    phen=$line
    echo $phen
    
    sumstat=${path}/sumstats_ldsc/${phen}_ldsc.sumstat
    out=${outdir}/${phen}_sumrhe

    python3 ${src}/sumrhe.py --pheno ${sumstat} \
                            --rhe ${demo}/rhe_trace_output_281k \
                            --all-snps \
                            --out ${out}

#done < "$path/phenos.txt"
#done < "$path/bili.txt"
done < "$path/residual_pheno_list.txt"


                  
                            #--trace ${demo}/trace_sum_281k_25.tr \


