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
parser.add_argument("--mem-eff", action='store_true', default=False,\
                    help='If --mem-eff flag is set, then we try to run in memory efficient mode (might be slightly slower)')


## TODO: implement memory-efficient flag (i.e., process one phenotype at a time)

class Sumrhe:
    def __init__(self, bim_path=None, rhe_path=None, sum_path=None, save_path=None, pheno_path=None, out=None, chisq_threshold=0, nbin=1, log=None, mem=False, allsnp=False):
        self.mem = mem
        self.log = log
        self.timezone = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo
        self.start_time = time.time()
        if (self.mem):
            self.log._log("Running with memory efficient mode")
            self.res = []
        self.log._log("Analysis started at: "+str(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self.start_time)))+" "+str(self.timezone))
        self.tr = trace.Trace(bimpath=bim_path, rhepath=rhe_path, sumpath=sum_path, savepath=save_path, log=self.log)
        if (rhe_path is not None):
            self.tr._read_all_rhe()
            if (save_path is not None):
                self.tr._save_trace()
        elif (sum_path is not None):
            self.tr._read_trace()
        self.snplist = self.tr.snplist
        self.nblks = self.tr.nblks
        if (allsnp):
            self.sums = sumstats.Sumstats(nblks=self.nblks, chisq_threshold=chisq_threshold, log=self.log)
        else:
            self.sums = sumstats.Sumstats(nblks=self.nblks, snplist = self.snplist, chisq_threshold=chisq_threshold, log=self.log)
        # check whether the path for pheno is a directory or a file
        if os.path.exists(pheno_path):
            if os.path.isdir(pheno_path):
                self.log._log("Reading phenotype sumstat files from a directory...")
                self.phen_dir = sorted([f for f in os.listdir(pheno_path) if f.endswith('.sumstat')])
                if (len(self.phen_dir) == 0):
                    self.log._log("!!! --pheno path has no valid phenotype summary files (.sumstat) !!!")
                    sys.exit(1)
                elif not self.mem:
                    for f in self.phen_dir:
                        self.sums._process(pheno_path.rstrip("/")+"/"+f)
            elif os.path.isfile(pheno_path):
                self.log._log("Reading a single phenotype sumstat file")
                self.sums._process(pheno_path)
            else:
                self.log._log("!!! --pheno path is invalid !!!")
                sys.exit(1)
        else:
            self.log._log("!!! --pheno path is invalid !!!")
            sys.exit(1)

        self.npheno = self.sums.npheno
        self.nsamp = self.sums.nsamp

        self.nbin = nbin # at the moment, SUMRHE only supports single-bin heritability
        self.nrows = nbin+1

        # the x, y arrays are reused each time
        self.xarray = np.zeros((self.nblks, self.nrows, self.nrows), dtype=float)
        self.yarray = np.zeros((self.nblks, self.nrows, 1), dtype=float)

        self.herits = np.zeros((self.npheno, self.nrows, self.nblks))
        self.hsums = np.zeros((self.npheno, self.nrows, 2))

        self.out = out

    def _solve_eqn(self, idx):
        for i in range(self.nblks):
            sol = np.linalg.solve(self.xarray[i], self.yarray[i])
            self.herits[idx][:, i] = sol.reshape((self.nrows, ))
        self.herits[idx] /= self.herits[idx].sum(axis=0)

    def _fill_normeq(self, idx):
        ''' fill up the xarray and yarray '''
        pred_tr = self.tr._calc_trace(self.nsamp[idx])
        for i in range(self.nblks):
            # TODO (for later): allow partitioned heritability
            self.xarray[i][0,0] = pred_tr[i]
            self.xarray[i][1,0] = self.nsamp[idx]
            self.xarray[i][0,1] = self.nsamp[idx]
            self.xarray[i][1,1] = self.nsamp[idx]

            self.yarray[i][0,0] = self.sums.yKys[idx][i]
            self.yarray[i][1,0] = self.nsamp[idx]

    def _run_jackknife(self, idx):
        ''' run snp-level block jackknife '''
        self.hsums[idx][:, 0] = self.herits[idx].mean(axis=1).reshape((self.nrows, ))
        self.hsums[idx][:, 1] = np.sqrt(np.sum(np.square(self.herits[idx] - self.hsums[idx][:,0].reshape((self.nrows, 1))), axis=1)*(self.nblks-1)/self.nblks)
        
    def _run(self):
        for i in range(self.npheno):
            self._fill_normeq(i)
            self._solve_eqn(i)
            self._run_jackknife(i)
        return self.hsums

    
    def _log(self):
        for i in range(self.npheno):
            h2 = self.hsums[i][0, 0]
            se = self.hsums[i][0, 1]
            self.log._log("^^^ Phenotype "+str(i)+" Estimated total heritability (h^2): "+format(h2, '.5f')+" SE: "+format(se, '.5f'))
        
        self.end_time = time.time()
        self.log._log("Analysis ended at: "+str(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(self.end_time)))+' '+str(self.timezone))
        self.log._log("run time: "+format(self.end_time - self.start_time, '.3f')+" s")
        if (self.out is not None):
            self.log._log("Saved log in "+ self.out + ".log")
            self.log._save_log(self.out+".log")
        return

if __name__ == '__main__':
    log = Logger()
    log._log(">>> SUM-RHE arguments")
    log._log("python3 sumrhe.py", end=" ")
    arg = sys.argv[1:]
    i = 0
    while i < len(arg):
        if arg[i].startswith('-'):
            if (i == 0):
                log._log(arg[i]+" "+arg[i + 1] if i + 1 < len(arg) else "")
            else:
                log._log("\t\t"+arg[i]+" "+arg[i + 1] if i + 1 < len(arg) else "")
            i += 1
        else:
            log._log(arg[i] if i==0 else '\t\t'+arg[i])
        i += 1
    args = parser.parse_args()
    if (args.trace is None) and (args.rhe is None):
        log._log("!!! Either trace sumamry or RHE trace output must be provided !!!")
        sys.exit(1)
    if (args.save_trace is not None) and (args.rhe is None):
        # TODO: allow combining the trace summaries with rhe traceoutputs
        log._log("!!! RHE trace output must be provided for --save-trace !!!")
    if (args.max_chisq is not None):
        if (args.max_chisq <= .0):
            log._log("!!! max-chisq must be a positive value !!!")
            sys.exit(1)

    sums = Sumrhe(bim_path=args.bim, rhe_path=args.rhe, sum_path=args.trace, save_path=args.save_trace, pheno_path=args.pheno,\
            chisq_threshold=args.max_chisq, log=log, out=args.out)
    sums._run()
    sums._log()

    


