# generate simulated phenotypes
split=10
h2=0.1
pol=_0.01

i=`expr $SGE_TASK_ID - 1`

dir=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${split}k

gen=${dir}/geno/${split}k_split
annot=${dir}/double_p0.05.txt
#annot=${dir}/singleannot.txt
src=/u/home/b/bronsonj/project-sriram/simulator/build_partition
par=${dir}/params/params_simul_${h2}${pol}${part}.txt
mafld=${dir}/maf_ld.txt
out=${dir}/sims/sims_${h2}_test${pol}${part}

mkdir -p ${out}

$src/Simulator -g $gen -annot $annot -simul_par $par -maf_ld $mafld -o $out/sim_${split}k_${i}  -k 10 -jn 100
