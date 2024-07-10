'''
Read the phenotype-specific summary statistics (generally PLINK)
The standard LDSC summary statistics format (.sumstat) works.
Format should be:
    SNPID, NMISS (or OBS_CT), Z
'''
import utils

import numpy as np
from os import listdir

class Sumstats:
    '''
    Sumstats class is called in the sums (sums.py) module. It's called after the trace module, which will then call this module
    in order to read in the phenotype sumstats. Current version can store only a single phenotype sumstats at a time (to save memory)
    '''
    def __init__(self, nblks=100, snplist = None, chisq_threshold=0, log=None, both_side=False, annot=None, nbins=1):
        self.log = log
        self.nblks = nblks
        self.nbins = nbins
        self.snplist = snplist
        if snplist is None:
            self.log._log("!!! Missing the list of SNPs used in trace calculation."+\
                  " All SNPs in the phenotype sumstats will be used. !!!")
        else:
            self.nsnps_trace = len(snplist)
        self.snpids = []
        self.annot = annot # TODO: convert the annot by matching with SNP IDs ("annot" passed through the arg is for trace SNPs)
        self.zscores = []
        self.zscores_bin = [] # list of list, where zscores are partitioned by bin assignment
        self.zscores_blk = [] # list of list of list, where zscores are partitioned by blk & bin assignment
        self.RHS = None # array of RHS for each pheno
        self.annot = annot # annotation for partitioned heritability (if None, assume single-bin)
        self.nsamp = 0
        self.nsnps = 0
        self.nsnps_blk = None # array for keeping track of number of SNPs in each leave-one-out blk
        self.nsnps_bin = None # total number of SNPs for each bin
        self.chisq_threshold = chisq_threshold # if positive value, then remove snps with chisq above it
        self.both_side = both_side # whether you filter SNPs on both sides. MUST have the snplist.
        self.matched_snps = None # array of matched SNPs
        if (self.both_side) and (snplist is None):
            self.log._log("!!! SNP list (.bim) must be input in order to perform both-side SNP filtering! !!!")
            return
        self.name = None
        self.removesnps = None
        
    def _filter_snps(self):
        ''' 
        TODO: Remove SNPs with chi-sq statistic above threshold. For 'problem' SNPs, change 
        beta to zero and reduce nsnps. This function must be called during match_snps
        return list of SNPs that were removed
        '''
        chisq = pow(self.zscores, 2)
        removesnps = []
        for i in range(self.nsnps):
            if (chisq[i] > self.chisq_threshold):
                self.zscores[i] = .0
                self.nsnps -= 1
                if (self.snplist is not None):
                    if (self.matched_snps[i] is not None):
                        removesnps.append(self.matched_snps[i])
        self.log._log("Removed "+str(len(removesnps))+" SNPs with chi-sq above the threshold")
        return removesnps
        
    def _read_sumstats(self, path, name):
        ''' Read in summary statistics for a single phenotype '''
        snpid = []
        Nmiss = []
        zscores = []
        self.name = name
        if (self.snplist is None):
            with open(path, 'r') as fd:
                next(fd)
                for line in fd:
                    val = line.split()
                    Nmiss.append(float(val[3]))
                    zscores.append(float(val[4]))
        else:
            with open(path, 'r') as fd:
                next(fd)
                for line in fd:
                    val = line.split()
                    snpid.append(str(val[0]))
                    Nmiss.append(float(val[3]))
                    zscores.append(float(val[4]))
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
        returns a nested list of (nblks, nbins, blk_size)
        '''
        matched_zscores = []
        nsnps_blk = np.zeros((self.nblks, self.nbins))
        if (self.snplist is None):
            # using all the SNPs
            for i in range(self.nblks):
                blk_size = self.nsnps//self.nblks
                blk_zscores = np.array(self.zscores[blk_size*i: blk_size*(i+1)] if (i < self.nblks-1) else self.zscores[blk_size*i: ])
                blk_annot = self.annot[blk_size*i:blk_size*(i+1)] if (i < self.nblks-1) else self.annot[blk_size*i:]
                partition, nsnps_partition = utils._partition_bin_non_overlapping(blk_zscores, blk_annot, self.nbins)
                matched_zscores.append(partition)
                nsnps_blk[i] = nsnps_partition
            self.log._log("Using " + str(self.nsnps) + " SNPs to calculate the RHS of the normal equation...")
        else:
            zscore_dict = dict(zip(self.snpids, self.zscores))
            snpid_dict = dict(zip(self.snpids, self.snplist))
            self.matched_snps = np.array([snpid_dict.get(pid) for pid in self.snpids])
            nmissing = len(self.snplist) - len(self.matched_snps)
            for i in range(self.nblks):
                blk_size = len(self.snplist)//self.nblks
                blk = self.snplist[blk_size*i: blk_size*(i+1)] if (i < self.nblks-1) else self.snplist[blk_size*i: ]
                blk_zscores = np.array([zscore_dict.get(pid, .0) for pid in blk])
                
                blk_annot = self.annot[blk_size*i:blk_size*(i+1)] if (i < self.nblks-1) else self.annot[blk_size*i:]
                partition, nsnps_partition = utils._partition_bin_non_overlapping(blk_zscores, blk_annot, self.nbins)
                matched_zscores.append(partition)
                nsnps_blk[i] = nsnps_partition
                # nmissing_blk = sum(1 for s in blk_zscores if s is None)
                # nmissing += nmissing_blk

            self.log._log("Matched "+str(len(self.snplist) - nmissing)+" SNPs in phenotype "+self.name+\
                  ", out of "+str(len(self.snplist))+" SNPs ("+str(nmissing)+" missing)")
        
        self.nsnps_blk = nsnps_blk
        self.zscores_blk = matched_zscores
        self.zscores_bin, self.nsnps_bin = utils._partition_bin_non_overlapping(self.zscores, self.annot, self.nbins)

        removesnps = None
        ## if SNP filtering is on
        if (self.chisq_threshold is not None):
            self.log._log("Filtering SNPs with chi-sq greater than "+str(self.chisq_threshold))
            self.removesnps = self._filter_snps()
            
        
    def _calc_RHS(self):
        ''' calculate the RHS of the normal equation '''
        self.rhs = np.full((self.nblks+1, self.nbins+1), self.nsamp)
        for i in range(self.nbins):
            total_zTz = np.dot(self.zscores_bin[i], self.zscores_bin[i])
            for j in range(self.nblks+1):
                if (j < self.nblks):
                    blk_zTz = np.dot(self.zscores_blk[j][i], self.zscores_blk[j][i])
                    self.rhs[j, i] = (total_zTz - blk_zTz)*self.nsamp/(self.nsnps_bin[i] - self.nsnps_blk[j][i])
                else:
                    self.rhs[j, i] = total_zTz*self.nsamp/self.nsnps_bin[i]
        self.log._log("Calculated the RHS for phenotype "+self.name)
        
    def _process(self, path, name):
        ''' 
        process a sumstat provided in the path
        '''
        self._read_sumstats(path, name)
        self._match_snps()
        self._calc_RHS()
        return self.removesnps