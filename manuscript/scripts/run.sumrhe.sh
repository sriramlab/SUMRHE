size=10k
tsize=281k
h2=0.1
#null
pol=_0.01

part=_p0.05

dir=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${size}
src=/u/home/b/bronsonj/project-sriram/sumrhe/src
out_dir=${dir}/sumrheout_${h2}_test${pol}${part}
trace=${dir}/ref${tsize}_sumtrace.tr

mkdir -p ${out_dir}
for i in {0..99}
do
    sumstat=${dir}/sums/sums_${h2}_test${pol}${part}_ldsc/sim_${size}_${i}.sumstat
    out=${out_dir}/${size}out_sim_${i}
    python3 ${src}/sumrhe.py --pheno $sumstat \
                      --out $out \
                      --trace $trace
done
