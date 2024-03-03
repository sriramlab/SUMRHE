'''
Read the phenotype-specific summary statistics (generally PLINK)
Format should be:
    SNPID, NMISS (or OBS_CT), Z
'''

import numpy as np
from os import listdir

class Sumstats:
    # Sumstats class is called in the main (sumrhe.py) module
    # it's called after the trace module, which will then call this module
    # in order to read in the phenotype sumstats.
    def __init__(self, nblks=100, snplist = None, chisq_threshold=0, log=None, both_side=False):
        self.log = log
        self.nblks = nblks
        self.snplist = snplist
        if snplist is None:
            self.log._log("!!! Missing the list of SNPs used in trace calculation."+\
                  " All SNPs in the phenotype sumstats will be used! !!!")
        else:
            self.nsnps_trace = len(snplist)
        self.snpids = []
        self.zscores = []
        self.zscores_blk = [] # list of list of arrays, where zscores are partitioned by blk assignment
        self.denom = [] # list of arrays, denominator for each pheno
        self.nsamp = []
        self.nsnps = []
        self.nsnps_blk = [] # list of arrays, keeping track of number of SNPs in each blk
        self.npheno = 0
        self.chisq_threshold = chisq_threshold # if positive value, then remove snps with chisq above it
        self.both_side = both_side # whether you filter SNPs on both sides. MUST have the snplist.
        self.matched_snps = [] # list of arrays
        if (self.both_side) and (snplist is None):
            self.log._log("!!! SNP list (.bim) must be input in order to perform both-side SNP filtering! !!!")
            return
        
    def _filter_snps(self, idx):
        ''' 
        TODO: Remove SNPs with chi-sq statistic above threshold. For 'problem' SNPs, change 
        beta to zero and reduce nsnps
        must be called during match_snps
        return list of SNPs that were removed
        '''
        chisq = pow(self.zscores[idx], 2)
        removesnps = []
        for i in range(self.nsnps[idx]):
            if (chisq[i] > self.chisq_threshold):
                self.zscores[idx][i] = .0
                self.nsnps[idx] -= 1
                if (self.matched_snps[idx][i] is not None):
                    removesnps.append(self.matched_snps[idx][i])
        return removesnps
        
    def _read_sumstats(self, path):
        ''' Read in summary statistics for a single phenotype '''
        snpid = []
        Nmiss = []
        zscores = []
        ## TEMPORARY FIX: Run with the old format (nmiss, betas, betas_se) w/o snplist
        ## TEMPORARY FIX: read in the LDSC sumstat
        if (self.snplist is None):
            with open(path, 'r') as fd:
                next(fd)
                for line in fd:
                    val = line.split()
                    Nmiss.append(float(val[3]))
                    #zscores.append(float(val[1])/float(val[2]))
                    zscores.append(float(val[4]))
        else:
            with open(path, 'r') as fd:
                next(fd)
                for line in fd:
                    val = line.split()
                    snpid.append(str(val[0]))
                    Nmiss.append(float(val[1]))
                    zscores.append(float(val[2]))
            self.snpids.append(snpid)
        nsamp = max(Nmiss)
        zscores = np.array(zscores)*np.sqrt(np.array(Nmiss)/nsamp)
        self.zscores.append(zscores)
        self.nsamp.append(nsamp)
        self.nsnps.append(len(zscores))
        self.npheno += 1
        return zscores
    
    def _match_snps(self, idx):
        ''' 
        match SNPs used for trace calculations & summary statistics
        returns a nested list of (nblks, blk_size)
        '''
        matched_zscores = []
        nsnps_blk = np.zeros(self.nblks)
        if (self.snplist is None):
            for i in range(self.nblks):
                blk_size = self.nsnps[idx]//self.nblks
                blk_zscores = np.array(self.zscores[idx][blk_size*i: blk_size*(i+1)] if (i < self.nblks-1) else self.zscores[idx][blk_size*i: ])
                matched_zscores.append(blk_zscores)
                nsnps_blk[i] = len(blk_zscores)
        else:
            zscore_dict = dict(zip(self.snpids[idx], self.zscores[idx]))
            snpid_dict = dict(zip(self.snpids[idx], self.snplist))
            matched_snps = np.array([snpid_dict.get(pid) for pid in self.snpids[idx]])
            self.matched_snps.append(matched_snps)
            nmissing = 0
            for i in range(self.nblks):
                blk_size = len(self.snplist)//self.nblks
                blk = self.snplist[blk_size*i: blk_size*(i+1)] if (i < self.nblks-1) else self.snplist[blk_size*i: ]
                blk_zscores = np.array([zscore_dict.get(pid) for pid in blk])
                matched_zscores.append(blk_zscores)
                nmissing_blk = sum(1 for s in blk_zscores if s is None)
                nmissing += nmissing_blk
                nsnps_blk[i] = len(blk_zscores) - nmissing_blk
            self.log._log("Matched "+str(len(self.snplist) - nmissing)+" SNPs in phenotype "+str(idx)+\
                  ", out of "+str(len(self.snplist))+" SNPs ("+str(nmissing)+" missing)")
        
        self.nsnps_blk.append(nsnps_blk)
        self.zscores_blk.append(matched_zscores)

        ## if SNP filtering is on
        if (self.chisq_threshold is not None):
            self._filter_snps(idx)
        return matched_zscores
            
        
    def _calc_denom(self, idx):
        ''' calculate denominator from the sumstats '''
        denom = np.zeros(self.nblks)
        total_zTz = np.dot(self.zscores[idx], self.zscores[idx])
        
        for i in range(self.nblks):
            blk_zscores = np.where(self.zscores_blk[idx][i] == None, 0, self.zscores_blk[idx][i])
            blk_zTz = np.dot(blk_zscores, blk_zscores)
            denom[i] = (total_zTz - blk_zTz)/(self.nsnps[idx] - self.nsnps_blk[idx][i]) - 1.0
        self.denom.append(denom)
        self.log._log("Calculated the denominator for phenotype "+str(idx))
        return denom
        
    def _process(self, path):
        ''' 
        process a sumstat provided in the path
        '''
        idx = self.npheno
        self._read_sumstats(path)
        self._match_snps(idx)
        self._calc_denom(idx)
