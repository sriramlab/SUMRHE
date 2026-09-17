# split genotype by CHR
. /u/local/Modules/default/init/modules.sh
module load plink

gen=/path/to/authorized_data/data/cal/filter4_no_mhc/filter4_no_mhc
out=/path/to/scratch/281k_chr
indv=/path/to/project/summaryRHE/300k/split_10k/281k.indv.txt

for chr in {1..22}
do
    plink --make-bed --bfile $gen  --chr ${chr}  --keep $indv --out ${out}/${chr}
done
