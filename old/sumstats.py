'''
Read the phenotype-specific summary statistics (generally PLINK)
Format should be:
    SNPID, NMISS (or OBS_CT), BETA, SE
'''

import numpy as np
from os import listdir

class Sumstats:
    # Sumstats class is called in the main (sumrhe.py) module
    # it's called after the trace module, which will then call this module
    # in order to read in the phenotype sumstats.
    def __init__(self, nblks=100, snplist = None, chisq_threshold=0, log=None):
        self.log = log
        self.nblks = nblks
        self.snplist = snplist
        if snplist is None:
            self.log._log("!!! Missing the list of SNPs used in trace calculation."+\
                  " All SNPs in the phenotype sumstats will be used! !!!")
        else:
            self.nsnps_trace = len(snplist)
        self.snpids = []
        self.betas = []
        self.betas_blk = [] # list of list of arrays, where betas are partitioned by blk assignment
        self.yKys = [] # list of arrays, yKy for each pheno
        self.nsamp = []
        self.nsnps = []
        self.nsnps_blk = [] # list of arrays, keeping track of number of SNPs in each blk
        self.npheno = 0
        self.chisq_threshold = chisq_threshold # if positive value, then remove snps with chisq above it
        
    def _filter_snps(self, threshold_chisq=None):
        ''' 
        TODO: Remove SNPs with chi-sq statistic above threshold. For 'problem' SNPs, change 
        beta to zero and reduce nsnps
        '''
        return
        
    def _read_betas(self, path):
        ''' Read in summary statistics for a single phenotype '''
        snpid = []
        Nmiss = []
        betas = []
        betas_se = []
        ## TEMPORARY FIX: Run with the old format (nmiss, betas, betas_se) w/o snplist
        if (self.snplist is None):
            with open(path, 'r') as fd:
                next(fd)
                for line in fd:
                    val = line.split()
                    Nmiss.append(float(val[0]))
                    betas.append(float(val[1]))
                    betas_se.append(float(val[2]))
        else:
            with open(path, 'r') as fd:
                next(fd)
                for line in fd:
                    val = line.split()
                    snpid.append(str(val[0]))
                    Nmiss.append(float(val[1]))
                    betas.append(float(val[2]))
                    betas_se.append(float(val[3]))
            self.snpids.append(snpid)
        nsamp = max(Nmiss)
        betas = np.array(betas)*np.array(Nmiss)/nsamp
        self.betas.append(betas)
        self.nsamp.append(nsamp)
        self.nsnps.append(len(betas))
        self.npheno += 1
        return betas
    
    def standardize_betas(self, path):
        ''' TODO: If the betas are not standardized - requires MAF of SNPs '''
        if (path is None) or (path == ""):
            self.log._log("!!! A valid MAF is required for standardization of the betas! !!!")
        ''' TODO: read maf '''
        for i in range(len(self.betas)):
            self.betas[i] *= np.sqrt(2*self.freq*(1-self.freq))
            self.betas_se[i] *= np.sqrt(2*self.freq*(1-self.freq))
            self.betas[i] /= np.sqrt(self.nsamp[i])*self.test.betas_se
            
    def _match_snps(self, idx):
        ''' 
        match SNPs used for trace calculations & summary statistics
        returns a nested list of (nblks, blk_size)
        must be called AFTER calling filter_SNPs (which will adjust the betas & the snp cnts)
        '''
        matched_betas = []
        nsnps_blk = np.zeros(self.nblks)
        if (self.snplist is None):
            for i in range(self.nblks):
                blk_size = self.nsnps[idx]//self.nblks
                blk_betas = np.array(self.betas[idx][blk_size*i: blk_size*(i+1)] if (i < self.nblks-1) \
                else self.betas[idx][blk_size*i: ])
                matched_betas.append(blk_betas)
                nsnps_blk[i] += len(blk_betas)
        else:
            beta_dict = dict(zip(self.snpids[idx], self.betas[idx]))
            nmissing = 0
            for i in range(self.nblks):
                blk_size = len(self.snplist)//self.nblks
                blk = self.snplist[blk_size*i: blk_size*(i+1)] if (i < self.nblks-1) \
                else self.snplist[blk_size*i: ]
                blk_betas = np.array([beta_dict.get(pid) for pid in blk])
                matched_betas.append(blk_betas)
                nmissing_blk = sum(1 for s in blk_betas if s is None)
                nmissing += nmissing_blk
                nsnps_blk[i] = len(blk_betas) - nmissing_blk
            self.log._log("Matched "+str(len(self.snplist) - nmissing)+" SNPs in phenotype "+str(idx)+\
                  ", out of "+str(len(self.snplist))+" SNPs ("+str(nmissing)+" missing)")
        
        self.nsnps_blk.append(nsnps_blk)
        self.betas_blk.append(matched_betas)
        return matched_betas
            
        
    def _calc_yKy(self, idx):
        ''' calculate yKy from the sumstats '''
        yKy = np.zeros(self.nblks)
        total_yXXy = np.dot(self.betas[idx], self.betas[idx])*pow(self.nsamp[idx], 2)
        
        for i in range(self.nblks):
            blk_betas = np.where(self.betas_blk[idx][i] == None, 0, self.betas_blk[idx][i])
            blk_yXXy = np.dot(blk_betas, blk_betas)*pow(self.nsamp[idx], 2)
            yKy[i] += (total_yXXy - blk_yXXy)/(self.nsnps[idx] - self.nsnps_blk[idx][i])
        self.yKys.append(yKy)
        self.log._log("Calculated yKy for phenotype "+str(idx))
        return yKy
        
    def _process(self, path):
        ''' 
        process a sumstat provided in the path
        '''
        idx = self.npheno
        self._read_betas(path)
        self._match_snps(idx)
        self._calc_yKy(idx)
