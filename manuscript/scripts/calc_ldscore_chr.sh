#!/bin/sh
source /u/local/Modules/default/init/modules.sh

chr=`expr $SGE_TASK_ID`

src_dir=/path/to/project/ldsc
all_bfile=/path/to/scratch/281k_chr
ldscore=/path/to/project/summaryRHE/300k/split_10k/ldscores_281k

# calculate ld score without annotation
python "$src_dir"/ldsc.py --l2 \
      --bfile $all_bfile/${chr} \
      --ld-wind-kb 2000 \
      --out $ldscore/${chr}
