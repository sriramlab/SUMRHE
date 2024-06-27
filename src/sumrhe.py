from logger import Logger
from gw_ldscore import GenomewideLDScore
from sums import Sumrhe

import argparse
import sys
import numpy as np
import utils

parser = argparse.ArgumentParser(description='SUM-RHE')
parser.add_argument("--trace", default=None, type=str, \
                    help='File path for trace summary statistics (.tr)')
parser.add_argument("--rhe", default=None, type=str, \
                   help='Directory where the rhe trace outputs (from same population & SNP sets) are stored.'
                    ' Each output must have .trace and .MN files')
parser.add_argument("--save-trace", default=None, type=str, \
                    help='File path for saving trace summaries calculated from RHE trace outputs (.tr)')
parser.add_argument("--pheno", default=None, type=str, \
                   help='File path for phenotype-specific summary statistics.'
                    ' If the path is a directory, all summary statistics (ending with .sumstat) will be used.')
parser.add_argument("--bim", default=None, type=str, \
                    help='File path for the reference .bim file used for trace calculation')
parser.add_argument("--max-chisq", action='store', default=None, type=float, \
                    help='Filter out SNPs with chi-sq statistic above the threshold.'
                    ' This can be done either only on the yKy or on both sides (use --filter-both-sides);'
                    ' with many non-polygenic SNPs, one-sided filtering might not be accurate')
parser.add_argument("--filter-both-sides", action='store_true', default=False,\
                    help='When filtering SNPs, remove their effects on both trace and yKy.'
                    ' This requires the (truncated) LD scores of all the SNPs used in trace calculation')
parser.add_argument("--ldscores", default=None, type=str, \
                    help='File path for LD scores of the reference SNPs. You may use either the traditional (truncated) LD scores (.l2.ldscore.gz) or genome-wide stochastic LD scores (.gw.ldscore.npy)')
parser.add_argument("--out", default=None, type=str, \
                    help='Output file path to save the analysis log and result (.log) or the genome-wide LD scores (.gw.ldscore.npy)')
parser.add_argument("--all-snps", action='store_true', default=False,\
                    help="Use all the SNPs in the phenotype sumamry statistics. Make sure this is safe to do so.")
parser.add_argument("--verbose", action="store_true", default=False,\
                    help='Verbose mode: print out the normal equations')
parser.add_argument("--suppress", action="store_true", default=False,\
                    help='Suppress mode: do not print out the outputs to stdout (log file only)')
parser.add_argument("--njack", default=100, type=int, \
                    help='Number of jackknife blocks (only if using ld projection matrix)')
parser.add_argument("--annot", default=None, type=str, \
                    help='Path of the annotation file (only if using partitioned heritability)')
parser.add_argument("--geno", default=None, type=str, \
                    help='Path of the genotype file to calculate the genome-wide LD scores. Calculates partitioned scores if --annot is also specified.')    

if __name__ == '__main__':
    args = parser.parse_args()
    log = Logger(suppress = args.suppress)
    log._log(">>> SUM-RHE arguments")
    log._log("python3 sumrhe.py", end=" ")
    arg = sys.argv[1:]
    i = 0
    while i < len(arg):
        if arg[i].startswith('-'):
            if (i == 0):
                log._log("\t"+arg[i]+" "+arg[i + 1] if i + 1 < len(arg) and not arg[i+1].startswith('-') else "")
                i += 1
            elif (i + 1 < len(arg)):
                if arg[i+1].startswith('-'):
                    log._log("\t"+arg[i])
                else:
                    log._log("\t"+arg[i]+" "+arg[i + 1])
                    i += 1
        else:
            log._log(arg[i] if i==0 else '\t\t'+arg[i])
        i += 1

    if (args.geno is not None):
        if (args.out is None):
            log._log("!!! An output path to save the genome-wide LD scores must be provided !!!")
            sys.exit(1)
        gwld = GenomewideLDScore(args.geno, args.annot, args.out, log)
        gwld._compute_ldscore()

    else:
        if (args.trace is None) and (args.rhe is None) and (args.ldscores is None):
            log._log("!!! Either trace sumamry, RHE trace output or LD score (truncated or stochastic) must be provided !!!")
            sys.exit(1)
        if (args.save_trace is not None) and ((args.rhe is None) and (args.ldscores is None)):
            # TODO: allow combining the trace summaries with rhe traceoutputs
            log._log("!!! RHE trace output or LD projection matrix must be provided for --save-trace !!!")
        if (args.max_chisq is not None):
            if (args.max_chisq <= .0):
                log._log("!!! max-chisq must be a positive value !!!")
                sys.exit(1)
        sums = Sumrhe(bim_path=args.bim, rhe_path=args.rhe, sum_path=args.trace, save_path=args.save_trace, pheno_path=args.pheno,\
            chisq_threshold=args.max_chisq, log=log, out=args.out, allsnp=args.all_snps, verbose=args.verbose, \
                filter_both=args.filter_both_sides, ldscores=args.ldscores, njack=args.njack, annot=args.annot)
        sums._run()
        sums._logoff()