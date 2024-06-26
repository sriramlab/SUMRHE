from sumstats import Sumstats
from trace import Trace
from logger import Logger
from gw_ldscore import GenomewideLDScore

import argparse
import os
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

class Sumrhe:
    def __init__(self, bim_path=None, rhe_path=None, sum_path=None, save_path=None, pheno_path=None, out=None, chisq_threshold=0, \
            log=None, mem=False, allsnp=False, verbose=False, filter_both=False, ldscores=None, njack=None, annot=None):
        self.mem = mem
        self.log = log
        self.start_time = utils._get_time()
        self.log._log("Analysis started at: "+utils._get_timestr(self.start_time))
        self.tr = Trace(bimpath=bim_path, rhepath=rhe_path, sumpath=sum_path, savepath=save_path, ldscores=ldscores, log=self.log, nblks=njack, annot=annot)
        self.snplist = self.tr.snplist
        self.nblks = self.tr.nblks
        self.annot = self.tr.annot
        self.nbins = self.tr.nbins
        if (allsnp):
            self.sums = sumstats.Sumstats(nblks=self.nblks, chisq_threshold=chisq_threshold, log=self.log, annot=self.annot, nbins=self.nbins)
        else:
            self.sums = sumstats.Sumstats(nblks=self.nblks, snplist = self.snplist, chisq_threshold=chisq_threshold, log=self.log, annot=self.annot, nbins=self.nbins)
        self.phen_dir = None
        # check whether the path for pheno is a directory or a file (or even regex). count # of phenotypes
        if os.path.exists(pheno_path):
            if os.path.isdir(pheno_path):
                self.log._log("Reading phenotype sumstat files from a directory...")
                phen_files = sorted([f for f in os.listdir(pheno_path) if f.endswith('.sumstat')])
                self.phen_dir = [pheno_path.rstrip("/")+"/"+name for name in phen_files]
                if (len(self.phen_dir) == 0):
                    self.log._log("!!! --pheno path has no valid phenotype summary files (.sumstat) !!!")
                    sys.exit(1)
            elif os.path.isfile(pheno_path):
                self.log._log("Reading a single phenotype sumstat file")
                self.phen_dir = [pheno_path]
            else:
                self.log._log("!!! --pheno path is invalid !!!")
                sys.exit(1)
        else:
            self.log._log("!!! --pheno path is invalid !!!")
            sys.exit(1)

        self.npheno = len(self.phen_dir)
        self.nsamp = []
        self.phen_names = [os.path.basename(name)[:-8] for name in self.phen_dir]

        self.herits = np.zeros((self.npheno, self.nblks+1, self.nbins+2)) # jackknife subsampled partitioned h2
        self.hsums = np.zeros((self.npheno, self.nbins+2, 2)) # partitioned h2 + total h2

        self.out = out
        self.filter_both = filter_both
        self.verbose = verbose

    def solve_linear_equation(self, X, y, method='lstsq'):
        '''
        Solve system of linear equations (either least square or QR)
        '''
        if (method == 'lstsq'):
            return np.linalg.lstsq(X, y, rcond=None)[0]
        else:
            Q, R = scipy.linalg.qr(X)
            return scipy.linalg.solve_triangular(R, np.dot(Q.T, y))

    def _calc_h2(self, idx):
        rhs = self.sums.rhs
        pred_tr = self.tr._calc_trace(self.nsamp[idx])
        for i in range(self.nblks+1):
            h2_est = self.solve_linear_equation(pred_tr[i], rhs[i])
            self.herits[idx][i] = np.append(h2_est, h2_est[:-1].sum()) # append total h2 as the last val
        if (self.verbose):
            sigmas = ["sigma^2_g" + str(i) for i in range(self.nbins)] + ["sigma^2_e"]
            sigmas_str = "["+", ".join(sigmas)+"]\n"
            self.log._log("Normal equation:\n"+np.array2string(pred_tr[self.nblks], precision=2, separator=', ')+"\n\t\t*\n"\
                +sigmas_str+"\t\t=\n"+np.array2string(rhs[self.nblks], precision=2, separator=', '))
            self.log._log("Solution:\n"+np.array2string(self.herits[idx][self.nblks][:-1], precision=3, separator=', '))
        return self.herits[idx]

    def _run_jackknife(self, idx):
        ''' run snp-level block jackknife '''
        self.hsums[idx, :, 0] = self.herits[idx, self.nblks] # h2 from all snps
        self.hsums[idx, :, 1] = np.sqrt(np.sum(np.square(self.herits[idx, :-1, :] - self.hsums[idx, :, 0].T), axis=0)*(self.nblks-1)/self.nblks) # jackknife SE

        
    def _run(self):
        for i in range(self.npheno):
            pheno_path=self.phen_dir[i]
            removesnps = self.sums._process(pheno_path, self.phen_names[i])
            self.tr._reset()
            if (removesnps is not None) and (self.filter_both):
                self.tr._filter_snps(removesnps)
            self.nsamp.append(self.sums.nsamp)
            self._calc_h2(i)
            self._run_jackknife(i)
        return self.hsums

    
    def _logoff(self):
        for i in range(self.npheno):
            h2 = self.hsums[i, -1][0]
            se = self.hsums[i, -1][1]
            self.log._log("^^^ Phenotype "+str(i)+" Estimated total heritability (h^2): "+format(h2, '.5f')+" SE: "+format(se, '.5f'))
        
        self.end_time = utils._get_time()
        self.log._log("Analysis ended at: "+utils._get_timestr(self.end_time))
        self.log._log("run time: "+format(self.end_time - self.start_time, '.3f')+" s")
        if (self.out is not None):
            self.log._log("Saved log in "+ self.out + ".log")
            self.log._save_log(self.out+".log")
        return

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