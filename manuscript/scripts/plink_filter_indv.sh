# create a subsample from genotype file
. /u/local/Modules/default/init/modules.sh
module load plink

gen=/u/project/sgss/UKBB/data/cal/filter4_no_mhc/filter4_no_mhc
indv=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_10k/281k.indv.txt
#indv=/u/home/b/bronsonj/project-sriram/summaryRHE/300k/split_10k/10k.indv.txt
out=/u/scratch/b/bronsonj/ref_geno/281k.10k_split

plink --make-bed --bfile $gen  --keep  $indv  --out $out
