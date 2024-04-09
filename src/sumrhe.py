import sumstats
import trace
from logger import Logger

import argparse
import os
import sys
import numpy as np
import time
import datetime

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
                    ' If the path is a directory, all summary statistics (ending with .sumstat) will be used.',
                    required=True)
parser.add_argument("--bim", default=None, type=str, \
                    help='File path for the reference .bim file used for trace calculation')
parser.add_argument("--max-chisq", action='store', default=None, type=float, \
                    help='Filter out SNPs with chi-sq statistic above the threshold.'
                    ' This can be done either only on the yKy or on both sides (use --filter-both-sides);'
                    ' with many non-polygenic SNPs, one-sided filtering might not be accurate')
parser.add_argument("--filter-both-sides", action='store_true', default=False,\
                    help='When filtering SNPs, remove their effects on both trace and yKy.'
                    ' This requires the (truncated) LD scores of all the SNPs used in trace calculation')
parser.add_argument("--ldscore", default=None, type=str, \
                    help='File path for LD scores of the reference SNPs, used for filtering non-polygenic SNPs')
parser.add_argument("--out", default=None, type=str, \
                    help='Output file path to save the analysis log and result (.log)')
parser.add_argument("--all-snps", action='store_true', default=False,\
                    help="Use all the SNPs in the phenotype sumamry statistics. Make sure this is safe to do so.")
parser.add_argument("--verbose", action="store_true", default=False,\
                    help='Verbose mode: print out the normal equations')
parser.add_argument("--suppress", action="store_true", default=False,\
                    help='Suppress mode: do not print out the outputs to stdout (log file only)')
parser.add_argument("--ldproj", default=None, type=str, \
                    help='File path for LD projection matrix in binary (.brz)')
parser.add_argument("--njack", default=100, type=int, \
                    help='Number of jackknife blocks (only if using ld projection matrix)')

class Sumrhe:
    def __init__(self, bim_path=None, rhe_path=None, sum_path=None, save_path=None, pheno_path=None, out=None, chisq_threshold=0, nbin=1, \
            log=None, mem=False, allsnp=False, verbose=False, filter_both=False, ldproj=None, njack=None):
        self.mem = mem
        self.log = log
        self.timezone = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo
        self.start_time = time.time()
        self.log._log("Analysis started at: "+str(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self.start_time)))+" "+str(self.timezone))
        self.tr = trace.Trace(bimpath=bim_path, rhepath=rhe_path, sumpath=sum_path, savepath=save_path, ldproj=ldproj, log=self.log, nblks=njack)
        if (rhe_path is not None):
            self.tr._read_all_rhe()
            if (save_path is not None):
                self.tr._save_trace()
        elif (sum_path is not None):
            self.tr._read_trace()
        elif (ldproj is not None):
            self.tr._read_ldproj()
        self.snplist = self.tr.snplist
        self.nblks = self.tr.nblks
        if (allsnp):
            self.sums = sumstats.Sumstats(nblks=self.nblks, chisq_threshold=chisq_threshold, log=self.log)
        else:
            self.sums = sumstats.Sumstats(nblks=self.nblks, snplist = self.snplist, chisq_threshold=chisq_threshold, log=self.log)
        self.phen_dir = None
        # check whether the path for pheno is a directory or a file. count # of phenotypes
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

        self.nbin = nbin # at the moment, SUMRHE only supports single-bin heritability
        self.nrows = nbin+1

        self.herits = np.zeros((self.npheno, self.nblks+1))
        self.hsums = np.zeros((self.npheno, self.nrows))

        self.out = out
        self.filter_both = filter_both
        self.verbose = verbose

    def _calc_h2(self, idx):
        denom = self.sums.denom
        pred_tr = self.tr._calc_trace(self.nsamp[idx])
        for i in range(self.nblks+1):
            numer = pred_tr[i]/self.nsamp[idx]-1.0
            self.herits[idx][i] = denom[i]/numer
        return self.herits[idx]

    def _run_jackknife(self, idx):
        ''' run snp-level block jackknife '''
        #self.hsums[idx][0] = self.herits[idx].mean()
        self.hsums[idx][0] = self.herits[idx, self.nblks]
        self.hsums[idx][1] = np.sqrt(np.sum(np.square(self.herits[idx][:-1] - self.hsums[idx][0]))*(self.nblks-1)/self.nblks)
        
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

    
    def _log(self):
        for i in range(self.npheno):
            h2 = self.hsums[i][0]
            se = self.hsums[i][1]
            self.log._log("^^^ Phenotype "+str(i)+" Estimated total heritability (h^2): "+format(h2, '.5f')+" SE: "+format(se, '.5f'))
        
        self.end_time = time.time()
        self.log._log("Analysis ended at: "+str(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self.end_time)))+' '+str(self.timezone))
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
                log._log("\t\t"+arg[i]+" "+arg[i + 1] if i + 1 < len(arg) else "")
            else:
                log._log("\t\t"+arg[i]+" "+arg[i + 1] if i + 1 < len(arg) else "")
            i += 1
        else:
            log._log(arg[i] if i==0 else '\t\t'+arg[i])
        i += 1
    if (args.trace is None) and (args.rhe is None) and (args.ldproj is None):
        log._log("!!! Either trace sumamry, RHE trace output or LD projection matrix must be provided !!!")
        sys.exit(1)
    if (args.save_trace is not None) and ((args.rhe is None) and (args.ldproj is None)):
        # TODO: allow combining the trace summaries with rhe traceoutputs
        log._log("!!! RHE trace output or LD projection matrix must be provided for --save-trace !!!")
    if (args.max_chisq is not None):
        if (args.max_chisq <= .0):
            log._log("!!! max-chisq must be a positive value !!!")
            sys.exit(1)

    sums = Sumrhe(bim_path=args.bim, rhe_path=args.rhe, sum_path=args.trace, save_path=args.save_trace, pheno_path=args.pheno,\
            chisq_threshold=args.max_chisq, log=log, out=args.out, allsnp=args.all_snps, verbose=args.verbose, \
                filter_both=args.filter_both_sides, ldproj=args.ldproj, njack=args.njack)
    sums._run()
    sums._log()