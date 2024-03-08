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
    # can contain only a single phenotype sumstats at a time (to save memory)
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
        self.zscores_blk = [] # list of arrays, where zscores are partitioned by blk assignment
        self.denom = None # array of denominator for each pheno
        self.nsamp = 0
        self.nsnps = 0
        self.nsnps_blk = None # array for keeping track of number of SNPs in each blk
        self.chisq_threshold = chisq_threshold # if positive value, then remove snps with chisq above it
        self.both_side = both_side # whether you filter SNPs on both sides. MUST have the snplist.
        self.matched_snps = None # array of matched SNPs
        if (self.both_side) and (snplist is None):
            self.log._log("!!! SNP list (.bim) must be input in order to perform both-side SNP filtering! !!!")
            return
        self.name = None
        
    def _filter_snps(self):
        ''' 
        TODO: Remove SNPs with chi-sq statistic above threshold. For 'problem' SNPs, change 
        beta to zero and reduce nsnps
        must be called during match_snps
        return list of SNPs that were removed
        '''
        chisq = pow(self.zscores, 2)
        removesnps = []
        for i in range(self.nsnps):
            if (chisq[i] > self.chisq_threshold):
                self.zscores[i] = .0
                self.nsnps -= 1
                if (self.matched_snps[i] is not None):
                    removesnps.append(self.matched_snps[i])
        return removesnps
        
    def _read_sumstats(self, path, name):
        ''' Read in summary statistics for a single phenotype '''
        snpid = []
        Nmiss = []
        zscores = []
        self.name = name
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
            self.snpids = snpid
        nsamp = max(Nmiss)
        zscores = np.array(zscores)*np.sqrt(np.array(Nmiss)/nsamp)
        self.zscores = zscores
        self.nsamp = nsamp
        self.nsnps = len(zscores)
        return zscores
    
    def _match_snps(self):
        ''' 
        match SNPs used for trace calculations & summary statistics
        returns a nested list of (nblks, blk_size)
        '''
        matched_zscores = []
        nsnps_blk = np.zeros(self.nblks)
        if (self.snplist is None):
            for i in range(self.nblks):
                blk_size = self.nsnps//self.nblks
                blk_zscores = np.array(self.zscores[blk_size*i: blk_size*(i+1)] if (i < self.nblks-1) else self.zscores[blk_size*i: ])
                matched_zscores.append(blk_zscores)
                nsnps_blk[i] = len(blk_zscores)
        else:
            zscore_dict = dict(zip(self.snpids, self.zscores))
            snpid_dict = dict(zip(self.snpids, self.snplist))
            matched_snps = np.array([snpid_dict.get(pid) for pid in self.snpids])
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
            self.log._log("Matched "+str(len(self.snplist) - nmissing)+" SNPs in phenotype "+self.name+\
                  ", out of "+str(len(self.snplist))+" SNPs ("+str(nmissing)+" missing)")
        
        self.nsnps_blk = nsnps_blk
        self.zscores_blk = matched_zscores

        removesnps = None
        ## if SNP filtering is on
        if (self.chisq_threshold is not None):
            removesnps = self._filter_snps()
        return removesnps
            
        
    def _calc_denom(self):
        ''' calculate denominator from the sumstats '''
        denom = np.zeros(self.nblks)
        total_zTz = np.dot(self.zscores, self.zscores)
        for i in range(self.nblks):
            blk_zscores = np.where(self.zscores_blk[i] == None, 0, self.zscores_blk[i])
            blk_zTz = np.dot(blk_zscores, blk_zscores)
            denom[i] = (total_zTz - blk_zTz)/(self.nsnps - self.nsnps_blk[i]) - 1.0
        self.denom = denom
        self.log._log("Calculated the denominator for phenotype "+self.name)
        return denom
        
    def _process(self, path, name):
        ''' 
        process a sumstat provided in the path
        '''
        self._read_sumstats(path, name)
        self._match_snps()
        self._calc_denom()