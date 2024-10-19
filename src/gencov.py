from sumstats import Sumstats
from trace import Trace
import utils

import numpy as np
import os
import sys

class GenCov:
    def __init__(self, bim_path=None, sum_path=None, save_path=None, pheno_path=None, out=None, chisq_threshold=0, \
            log=None, mem=False, allsnp=False, verbose=False, filter_both=False, ldscores=None, njack=None, annot=None):
        self.mem = mem
        self.log = log
        self.start_time = utils._get_time()
        self.log._log("Analysis started at: "+utils._get_timestr(self.start_time))
        self.tr = Trace(bimpath=bim_path, sumpath=sum_path, savepath=save_path, ldscores=ldscores, log=self.log, nblks=njack, annot=annot, verbose=verbose)
        self.snplist = self.tr.snplist
        self.nblks = self.tr.nblks
        self.annot = self.tr.annot
        self.annot_header = self.tr.annot_header
        self.nbins = self.tr.nbins
        if (allsnp):
            self.sums = Sumstats(nblks=self.nblks, chisq_threshold=chisq_threshold, log=self.log, annot=self.annot, nbins=self.nbins)
        else:
            self.sums = Sumstats(nblks=self.nblks, snplist = self.snplist, chisq_threshold=chisq_threshold, log=self.log, annot=self.annot, nbins=self.nbins)
        self.phen_dir = None
        # check whether the path for pheno is a directory or a file (or even regex). count # of phenotypes
        # TODO: allow regex matching for file names
        if os.path.exists(pheno_path):
            if os.path.isdir(pheno_path):
                self.log._log("Reading phenotype sumstat files from a directory...")
                phen_files = sorted([f for f in os.listdir(pheno_path) if f.endswith('.sumstat')])
                self.phen_dir = [pheno_path.rstrip("/")+"/"+name for name in phen_files]
                if (len(self.phen_dir) == 0):
                    self.log._log(f"!!! --pheno path {pheno_path} has no valid phenotype summary files (.sumstat) !!!")
                    sys.exit(1)
            elif os.path.isfile(pheno_path):
                self.log._log("Reading a single phenotype sumstat file")
                self.phen_dir = [pheno_path]
            else:
                self.log._log(f"!!! --pheno path {pheno_path} is invalid !!!")
                sys.exit(1)
        else:
            self.log._log(f"!!! --pheno path {pheno_path} is invalid !!!")
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
        if (self.verbose):
            self.log._log("Solution & jackknife SE:\n"+np.array2string(self.hsums, precision=5, separator=', '))

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
            if (self.nbins > 1):
                for j in range(self.nbins):
                    h2 = self.hsums[i, j, 0]
                    se = self.hsums[i, j, 1]
                    self.log._log(f"^^^ Phenotype {i} Estimated partitioned heritability bin {self.annot_header[j]} (h^2_{j}): {h2:.5f} SE: {se:.5f}")
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
