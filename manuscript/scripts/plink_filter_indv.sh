# create a subsample from genotype file
. /u/local/Modules/default/init/modules.sh
module load plink

gen=/path/to/authorized_data/data/cal/filter4_no_mhc/filter4_no_mhc
indv=/path/to/project/summaryRHE/300k/split_10k/281k.indv.txt
#indv=/path/to/project/summaryRHE/300k/split_10k/10k.indv.txt
out=/path/to/scratch/ref_geno/281k.10k_split

plink --make-bed --bfile $gen  --keep  $indv  --out $out
