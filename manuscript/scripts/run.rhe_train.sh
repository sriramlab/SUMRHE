# script for generating trace summaries

size=50k

split=10k
dir=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_${split}

gen=/u/scratch/b/bronsonj/ref_geno/281k.subset.${size}_split
src=/u/home/b/bronsonj/project-sriram/RHE-mc/build
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


