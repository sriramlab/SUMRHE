# Supplemental code

## Here are the main code and scripts used in this study.

## The ```simulator``` is a modified version of the lab phenotype simulator on Github (https://github.com/sriramlab/Simulator) that supports partitioned effect sizes & few other features.

## Here's a quick overview of the scripts:
### On simulations
1. Creating genotype
   
    ```plink_filter_indv.sh```: creates subset of the original ```filter4_no_mhc.bed``` file for the experiments.
   
    ```plink_split_chr.sh```: splits a large genotype (> 200k) by chromosomes, so that LD scores or SNP taggings can be calculated.

2. Simulating phenotype
   
    ```makesim.sh```: uses the ```simulator``` code to generate phenotypes of diverse genetic architectures.

3. Preparing input

    ```run.rhe_train.sh```: run ```RHE-mc``` on the reference genotype to generate trace summaries (using the ```-tr``` option)
   
    ```run.sumher_tagging.sh```: calculate ```SumHer``` SNP-taggings on the reference genotype
   
    ```calc_ldscore_chr.sh```: calculate LDSC LD scores for each chromosome of the reference genotype.
   
    ```gen_sumstats.sh```: run GWAS on the simulated phenotypes using ```PLINK 2.0```.

4. Running methods
   
    ```run.bolt_reml.sh```: run ```BOLT-REML``` to estimate heritability from genotype and simulated phenotypes.
   
    ```run.gcta_reml.sh```: run ```GCTA-REML``` to estimate heritability from genotype and simulated phenotypes.
   
    ```run.rhe.sh```: run ```RHE-mc``` to estimate heritability from genotype and simulated phenotypes.
   
    ```run.sumrhe.sh```: run ```SUM-RHE``` to estimate heritability from trace summaries (.tr) and PLINK GWAS summary statistics.
   
    ```run.ldsc.sh```: run ```LDSC``` to estimate heritability from (LDSC) LD scores and PLINK GWAS summary statistics.
   
    ```run.sumher.sh```: run ```SumHer``` to estimate heritability from SNP taggings and PLINK GWAS summary statistics.
   
### On real UKBB phenotypes
1. Preparing input
   
    ```regress_on_cov.ipynb```: covariates of top 20 PC's, sex and age were regressed out from the real phenotypes.
   
    ```residual_pheno_list.txt```: list of real phenotypes that were plotted in the paper (selected by the highest z-scores of heritability estimates)
   
    For summary-based methods, same scripts (but with in-sample genotype) were used to calculate input summaries (trace, LD score, SNP-taggings)
   
2. Running methods
   
    ```run.rhe.realphen.sh```, ```run.sumrhe.realphen.sh```, ```run.ldsc.realphen.sh```, ```run.sumher.realphen.sh```, ```run.sumher_ldak.realphen.sh```

### Additional notes:
Given ```RHE-mc```  is no longer maintained, we recommend using ```pyRHE``` or ```GENIE``` from our lab, which has all the functionalities.

