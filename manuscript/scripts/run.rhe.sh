split=10k
h2=0.1
pol=_0.01
dir=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${split}

part=
#_p0.05

gen=${dir}/geno/${split}_split
src=/u/home/b/bronsonj/project-sriram/RHE-mc/build_sums
out=${dir}/rheout_${h2}_test${pol}${part}
phen=${dir}/sims/sims_${h2}_test${pol}${part}
annot=${dir}/singleannot.txt

i=`expr $SGE_TASK_ID - 1`

mkdir -p $out

outfile=$out/singleout_${split}_run_${i}.txt


ls $outfile > /dev/null 2>&1

if [ $? -ne 0 ]
then
    #{ time $src/RHEmc_mem -g $gen -p $phen/sim_${split}k_${i}.phen -k 100 -jn 1000  -o $out/singleout_${split}k_run_${i}.txt -annot $annot; } 2> timing_test_${i}.txt
    $src/RHEmc_mem -g $gen -p $phen/sim_${split}_${i}.phen -k 100 -jn 1000  -o $out/singleout_${split}_run_${i}.txt -annot $annot
fi


