

source /u/local/Modules/default/init/modules.sh
module load anaconda2/2019.10


size=10k
h2=0.1
#null
pol=_0.01

part=_p0.05

dir=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${size}

src_dir=/u/home/b/bronsonj/project-sriram/ldsc
out_dir=${dir}/ldscout_${h2}_test${pol}${part}_refld
ld_dir=${dir}/ldscores_281k
 
mkdir -p ${out_dir}

for i in {0..99}
do
    sumstat=${dir}/sums/sums_${h2}_test${pol}${part}_ldsc/sim_${size}_${i}.sumstat
    out=${out_dir}/${size}out_sim_${i}

    python "$src_dir"/ldsc.py \
          --h2 $sumstat\
          --ref-ld ${ld_dir}/all \
          --w-ld ${ld_dir}/all \
          --out $out \
          --not-M-5-50 \
          --n-blocks 1000 \
          --chisq-max 99999
done
