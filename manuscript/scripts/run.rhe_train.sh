# script for generating trace summaries

size=50k

split=10k
dir=/path/to/project/summaryRHE/300k/split_${split}

gen=/path/to/scratch/ref_geno/281k.subset.${size}_split
src=/path/to/project/RHE-mc/build
out=${dir}/rheout_train_${size}
annot=${dir}/singleannot.txt

i=`expr $SGE_TASK_ID - 1`

mkdir -p $out

outfile=$out/singleout_${size}_run_${i}.txt


ls $outfile > /dev/null 2>&1

if [ $? -ne 0 ]
then
    $src/RHEmc_mem -g $gen -k 100 -jn 1000  -o $out/singleout_${size}_run_${i}.txt -annot $annot -tr
fi


