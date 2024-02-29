

size=50k
h2=0.1
pol=_0.01
#_0.01
#_0.01
#_0.1
pc=_p0.05
#_pstrat

if [ "$size" == "50k" ]
then
    tsize=241k
else
    tsize=281k
fi


phen_dir=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${size}
#phen=${phen_dir}/sums_0.25_241k_ldsc
#phen=${phen_dir}/sums/sums_${h2}_test${pol}${pc}_ldsc
phen=${phen_dir}/sums_${h2}_test${pol}${pc}_ldsc
dir=/u/home/b/bronsonj/project-sriram/sumrhe
src=${dir}/src
demo=${dir}/data
outdir=${phen_dir}/sumrheouts/sumsout_${h2}_test${pol}${pc}
#outdir=${phen_dir}/sumrheouts_train
#outdir=${phen_dir}/sumrheouts/sumsout_${h2}_test${pol}${pc}
#outdir=${phen_dir}/sumrheouts/ntrains/sumsout_${h2}_test${pol}${pc}_${tsize}_${n}

mkdir -p ${outdir}


#for i in {0..0};
for i in {0..99};
do
    python3 ${src}/sumrhe.py --pheno ${phen}/sim_${size}_${i}.sumstat \
                            --all-snps \
                            --trace ${demo}/trace_sum_281k_25.tr \
                            --out ${outdir}/sim_${size}_${i} 
    #> /dev/null
done
                  
                            #--rhe ${demo}/rhe_trace_output_241k \
                            #--save-trace ${demo}/trace_sum_241k_25 \


