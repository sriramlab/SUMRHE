# SUM-RHE
**SUM**mary statistics-based **R**andomized **H**aseman-**E**lston regression

## Prerequisites
Python (3.0 or higher) is required to run SUM-RHE.

## How to install/run:
Installing and running SUM-RHE is very simple. Clone the codebase and simply run with appropriate input/arguments.
```
git clone https://github.com/bronsonj98/SUMRHE.git
```

An example run command is as follows (run.sh in the ```example``` directory):
```
python3 ../src/sumrhe.py --pheno ./sim_50k_h_0.25_p_0.01.sumstat \
                  --ldscores ./double_uniform_0.2_ldscores.trunc.npy \
                  --annot ./double_uniform_0.2.txt \
                  --verbose \
                  --njack 1000
```

SUM-RHE support partitioned heritability as well. Currently, the annotation file used for calculating the ld scores must also be provided separately (to be updated).

If you'd like to create your own trace summaries, please refer to the PyRHE program from our lab: https://github.com/sriramlab/PyRHE.

You would need your own individual-level genotype for this (running with ```-tr``` option will print the trace summaries).


## Parameters

```
--trace : The path of aggregate trace summaries (.tr)
--rhe : The path of RHE trace summaries (.trace) and corresponding metadata (.MN)
--pheno : The path (either directory or file) with GWAS summary statistics (.sumstat)
--bim : File path for the reference .bim file used for trace calculation (optional)
--out : Output file path to save the analysis log and result (.log)
--save-trace : File path for saving trace summaries calculated from RHE trace outputs (.tr/.MN)
--max-chisq : Filter out SNPs with chi-sq statistic above the threshold.
--filter-both-sides : When filtering SNPs, remove their effects on both trace and yKy.
--ldscore : File path for LD scores of the reference SNPs. You may use either the traditional (truncated) LD scores (.l2.ldscore.gz) or genome-wide stochastic LD scores (.gw.ldscore.npy)
--all-snps : Use all the SNPs in the phenotype sumamry statistics. Make sure this is safe to do so.
--verbose : Verbose mode: print out the normal equations
--suppress : Suppress mode: do not print out the outputs to stdout (log file only)
--njack : Number of jackknife blocks (if using ld projection matrix)
--annot : Path of the annotation file (if using partitioned heritability)
```

## References
You may refer to the following bioRxiv preprint for more details on how SUM-RHE works & benchmark results:
https://www.biorxiv.org/content/10.1101/2024.03.09.584258v1.full
