
dir=/u/home/b/bronsonj/project-sriram/sumrhe
src=${dir}/src
demo=${dir}/data

python3 ${src}/sumrhe.py --pheno ${demo}/sumstats/10ksums_0.25_1.0/sim_0.sumstat \
                  --bim ${dir}/10k_split.bim \
                  --rhe ${demo}/rhe_trace_output_281k \
                  --save-trace ${demo}/sumtrace/trace_281k_ref \
                  --out ./out
                  


