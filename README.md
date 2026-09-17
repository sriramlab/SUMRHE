# SUM-RHE

**SUM**mary statistics-based **R**andomized **H**aseman-**E**lston regression

```sumrhe``` is an efficient tool for accurately estimating:
1. heritability of phenotypes from summary statistics
2. genome-wide LD scores from biobank data

## Data and examples

This repository distributes code only. Bundled research and example datasets have
been removed. The commands below are templates for inputs you are authorized to
use; replace `/path/to/...` with paths outside this checkout.

Do not commit participant identifiers, genotypes, phenotypes, covariates, or
research data and outputs. Simulating phenotypes does not make the underlying
genotypes synthetic. See [CONTRIBUTING.md](CONTRIBUTING.md) before adding files or
updating an older clone.

## How to get started
You may set up ```sumrhe``` using Conda (Anaconda3 or Miniconda) or virtual environment and pip (Miniconda or Python3)

### Using Conda (Anaconda3 or Miniconda)
1. Clone the repository
```bash
git clone https://github.com/sriramlab/SUMRHE.git
cd SUMRHE
```

2. Create the conda environment
```bash
conda env create -f environment.yml
```

3. Activate the conda environment
```bash
conda activate sumrhe
```

### Using Virtual Environment and Pip (Miniconda or Python3)
1. Clone the repository
```bash
git clone https://github.com/sriramlab/SUMRHE.git
cd SUMRHE
```

2. Create a virtual environment
```bash
python -m venv sumrhe
```

3. Activate the virtual environment
```bash
# on Windows
sumrhe\Scripts\activate
# on macOS/Linux
source sumrhe/bin/activate
```

4. Install dependencies
```bash
pip install -r requirements.txt
```

## How to use ```sumrhe```
### 1. Estimating partitioned heritability from summary statistics

```sumrhe``` can accurately (i.e., comparable to methods that use individual-level data) estimate heritability from summary-level data.
To estimate (partitioned) heritability, you need the following: 
1. LD scores (either genome-wide or fixed-window) **OR** trace summaries calculated from in-sample or reference population genotype.
2. GWAS summary statistics for the trait of interest

From the repository root, use your own summary statistics, LD scores, and annotation:
```
python3 src/sumrhe.py --pheno /path/to/inputs/trait.sumstat \
                  --ldscores /path/to/inputs/reference.gw.ldscore.gz \
                  --annot /path/to/inputs/annotation.txt \
                  --out /path/to/outputs/heritability \
                  --verbose \
                  --njack 1000
```
Currently, to estimate partitioned heritability, the "thin" annotation file used for calculating the LD scores must also be provided separately (to be updated).

Similarly, you may estimate heritability using trace summaries. Supply the prefix of your matching `.tr` and `.MN` files:
```
python3 src/sumrhe.py --pheno /path/to/inputs/trait.sumstat \
                  --trace /path/to/inputs/reference \
                  --annot /path/to/inputs/annotation.txt \
                  --out /path/to/outputs/heritability_trace \
                  --verbose
```

If you'd like to create your own trace summaries, please refer to the ```PyRHE``` program from our lab: https://github.com/sriramlab/PyRHE. You may also use the C++ version of ```GENIE```.
You would need your own individual-level genotype for this (running with ```-tr``` option will save the trace summaries). While less flexible than using SNP-level LD scores, the trace summaries are directly estimated with ```GENIE``` or ```pyRHE```, and the files are a lot smaller in size (less than 0.1 MB).

### 2. Estimating genome-wide LD scores

Other methods use sliding fixed-sized windows (typically < 2Mb) to estimate the LD scores of the SNPs. This results in under-estimation of LD scores, as the long-range LD (> 2Mb) is not captured. Often times, this results in upward-biased estimates of heritability. One of the advantages of Randomized Haseman-Elston regression is that it can efficiently estimate genome-wide correlations between SNPs through random projection. In the original RHE papers, this is used to estimate the trace of squared kinship matrix. ```sumrhe``` extends this idea further by estimating the genome-wide (partitioned) LD scores.

Use your own PLINK `.bed`, `.bim`, and `.fam` files, passing their common prefix:
```
python3 src/sumrhe.py --geno /path/to/inputs/genotypes \
                  --annot /path/to/inputs/annotation.txt \
                  --out /path/to/outputs/reference \
                  --nvecs 100 \
                  --nworkers 8
```
This command creates `reference.gw.ldscore.gz` and `reference.log` in your output directory. Runtime depends on the input size. The first three LD-score columns are metadata ('CHR', 'SNP', 'BP'), and the remaining columns contain the (partitioned) LD scores.

## Parameters

```
--trace : File path for trace summary statistics (.tr) and corresponding metadata (.MN). If the path is a directory, all trace summaires (ending with .tr) will be used by aggregating them.
--save-trace : File path for saving (aggregated) trace summaries (.tr) and corresponding metadata (.MN)
--pheno : File path for phenotype-specific summary statistics (.sumstat). If the path is a directory, all summary statistics (ending with .sumstat) will be used.
--bim : File path for the reference .bim file used for trace calculation (optional)
--out : Output file path to save the analysis log and result (.log) or the genome-wide LD scores (.gw.ldscore.gz)
--max-chisq : Filter out SNPs with chi-sq statistic above the threshold.
--filter-both-sides : When filtering SNPs, remove their effects on both trace and yKy.
--ldscore : File path for LD scores of the reference SNPs. You may use either the traditional (truncated) LD scores (.l2.ldscore.gz) or genome-wide stochastic LD scores (.gw.ldscore.gz)
--all-snps : Use all the SNPs in the phenotype sumamry statistics. Make sure this is safe to do so.
--verbose : Verbose mode: print out the normal equations
--suppress : Suppress mode: do not print out the outputs to stdout (log file only)
--njack : Number of jackknife blocks (only if using LD scores as input)
--annot : Path of the annotation file (if using partitioned heritability)
--geno : Path of the genotype file to calculate the genome-wide LD scores. Calculates partitioned scores if --annot is also specified.
--nworkers : Number of workers for multiprocessing to calculate stochastic genome-wide LD scores. Default is 4.
--nvecs : Number of random vectors to use for estimating stochastic genome-wide LD scores. Default is 10.
--step_size : Number of SNPs to process in each step of estimating stochastic genome-wide LD scores. Default is 1000.
--seed : Seed for estimating stochastic genome-wide LD scores. If not specified, the default numpy (pseudo) random number generator will be used.
```

## TODO's
✅ partitioned heritability

✅ stochastic genome-wide LD scores

✅ better SE estimates with LD scores

✅ easier input file formatting

☑️ both-side filtering of outlier SNPs

☑️ genetic correlation

## References
```SUM-RHE``` is now published on Genome Research as an open access article. For references, please cite
```
Jeong, M., Pazokitoroudi, A., Liu, Z., & Sankararaman, S. (2024). Scalable summary statistics-based heritability estimation method with individual genotype level accuracy. Genome research, gr.279207.124. Advance online publication. https://doi.org/10.1101/gr.279207.124
```

