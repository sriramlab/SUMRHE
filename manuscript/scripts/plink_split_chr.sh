# split genotype by CHR
. /u/local/Modules/default/init/modules.sh
module load plink

gen=/u/project/sgss/UKBB/data/cal/filter4_no_mhc/filter4_no_mhc
out=/u/scratch/b/bronsonj/281k_chr
indv=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_10k/281k.indv.txt

for chr in {1..22}
do
    plink --make-bed --bfile $gen  --chr ${chr}  --keep $indv --out ${out}/${chr}
done
