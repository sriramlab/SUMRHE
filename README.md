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
python3 sumrhe.py --pheno ./test --all-snps --trace ${demo}/trace_sum_281k_25.tr --max-chisq 15
```

## Parameters

```
--trace : The path of aggregate trace summaries (.tr)
--rhe : The path of RHE trace summaries (.trace) and corresponding metadata (.MN)
--pheno : The path (either directory or file) with GWAS summary statistics (.sumstat)
--bim : File path for the reference .bim file used for trace calculation
--out : Output file path to save the analysis log and result (.log)
--save-trace : File path for saving trace summaries calculated from RHE trace outputs (.tr)
--max-chisq : Filter out SNPs with chi-sq statistic above the threshold.
--filter-both-sides : When filtering SNPs, remove their effects on both trace and yKy.
--ldscore : File path for LD scores of the reference SNPs, used for filtering non-polygenic SNPs.
--all-snps : Use all the SNPs in the phenotype sumamry statistics. Make sure this is safe to do so.
--verbose : Verbose mode: print out the normal equations
--suppress : Suppress mode: do not print out the outputs to stdout (log file only)
--ld_proj : File path for LD projection matrix in binary (.brz)
--njack : Number of jackknife blocks (only if using ld projection matrix)
```

## References
You may refer to the following bioRxiv preprint for more details on how SUM-RHE works & benchmark results:
https://www.biorxiv.org/content/10.1101/2024.03.09.584258v1.full
