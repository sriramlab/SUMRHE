'''
Read in the trace (or LD) summary statistics to estimate the trace for the summary stats.
For .tr file, there should be (njackknife+1) rows excluding the header, where
the last row is the trace summary (i.e., sum of the LD scores) from the entire genotype.
The SNP sets should be the same for all the outputs -- it is recommended to input .bim for the list of SNPs (read in once)

For .ldscore.gz, standard LDSC ld score format works.
'''
import utils

import numpy as np
import pandas as pd
from os import listdir, path
import sys

class Trace:
    def __init__(self, bimpath = None, sumpath=None, savepath=None, log=None, ldscores=None, nblks=100, annot=None, verbose=False):
        self.log = log
        self.sumpath = sumpath
        self.savepath = savepath
        self.ldscorespath = ldscores
        self.ldscores = None
        self.sums = []
        self.nblks = nblks # nblks specified only if using ld proj; otherwise it'll be overwritten by trace summaries.
        self.ntrace = 0
        self.K = []
        self.snplist = None
        self.nsamp = [] # number of samples used for trace summaries; can be of varying size
        self.nsnps = 0 # number of SNPs is fixed
        self.nsnps_blk = None # array for keeping track of number of SNPs in each leave-one-out blk
        self.nsnps_bin = None
        self.nsnps_blk_filt = None # this is the nsnps for each jn, which is filtered in case of --filter-both-sides
        self.verbose = verbose
        if (bimpath is None) or (bimpath == ""):
            if (self.ldscorespath is None):
                self.log._log("!!! SNP list (.bim) is generally recommended if using trace summaries (.tr) !!!")
        elif (bimpath.endswith(".bim")):
            with open(bimpath, 'r') as fd:
                for line in fd:
                    self.snplist.append(line.split()[1])
        else:
            self.log._log(f'!!! {bimpath} is not a .bim file! !!!')
        
        if (sumpath is not None):
            self._read_all_trace()
            if (savepath is not None):
                self._save_trace()
        elif (ldscores is not None):
            self._read_ldscores()
        
        self._read_annot(annot)
    
    def _read_annot(self, annot_path):
        if (annot_path is None):
            self.annot_header = np.array(['L2'])
            self.annot = np.ones((self.nsnps, 1))
            self.log._log("Running with single component")
        else:
            self.annot_header, self.annot = utils._read_with_optional_header(annot_path)
            if (self.annot.ndim == 1):
                self.annot = self.annot.reshape(-1, 1)
            if (self.annot_header is None):
                self.annot_header = np.array(['L2_'+str(i) for i in range(self.annot.shape[1])])
            self.log._log("Read SNP partition annotation of dimensions "+str(self.annot.shape))
        if (self.nbins != self.annot.shape[1]) or (self.nsnps != self.annot.shape[0]):
            self.log._log("!!! number of components in annotation does not match the input trace summary !!!")
            sys.exit(1)
        self.blk_size = self.nsnps//self.nblks
        self.log._log("Number of jackknife blocks: "+str(self.nblks)+", blk_size: "+str(self.blk_size))
        self.nsnps_bin = self.annot.sum(axis=0)
        self.nsnps_blk = np.full((self.nblks+1, self.nbins), self.nsnps_bin)
        for i in range(self.nblks):
            idx_start = self.blk_size*i
            idx_end = self.blk_size*(i+1)
            if (i==self.nblks-1):
                idx_end = self.nsnps
            self.nsnps_blk[i] -= self.annot[idx_start: idx_end].sum(axis=0)
       
    
    def _save_trace(self):
        ''' Save trace summaries as a file'''
        with open(self.savepath+".MN", 'w') as fd:
            fd.write("NSAMPLE,NSNPS,NBLKS,NBINS,K\n")
            fd.write(f"{self.nsamp:.0f},{self.nsnps:.0f},{self.nblks:.0f},{self.nbins:.0f},{self.effective_K:.0f}")

        with open(self.savepath+".tr", 'w') as fd:
            header_str = ','.join(f'LD_SUM_{i:d}' for i in range(self.nbins))
            fd.write(header_str+",NSNPS_JACKKNIFE\n")
            for j in range(self.nblks+1):
                for k in range(self.nbins):
                    row_str = ','.join(f'{self.sums[j,k,l]:.3f}' for l in range(self.nbins))
                    row_str += f',{self.nsnps_blk[j, k]:.0f}\n'
                    fd.write(row_str)
        self.log._log(f"Saved trace summary into {self.savepath}(.tr/.MN)")
    
    def _read_trace(self, filename, idx):
        '''
        Read trace summaries (block-wise LD scores)
        '''
        # read in metadata
        nsamp, nsnps, nblks, nbins = 0, 0, 0, 0
        with open(filename+".MN", 'r') as fd:
            next(fd)
            nsamp, nsnps, nblks, nbins, K = map(int, fd.readline().split(','))
        self.nsamp.append(nsamp)
        self.K.append(K)

        if not idx:
            self.nsnps = nsnps
            self.nblks = nblks
            self.nbins = nbins
        elif (nblks != self.nblks):
            self.log._log("!!! Trace summary "+filename+" has incorrect number of jackknife blocks !!!")
            return
        elif (nbins != self.nbins):
            self.log._log("!!! Trace summary "+filename+" has incorrect number of annotation bins !!!")
            return

        sums = np.zeros((self.nblks+1, self.nbins, self.nbins))
        nsnps_blk = np.zeros((self.nblks+1, self.nbins))
        
        # read in trace values
        for cnt, vals in enumerate(utils._read_multiple_lines(filename+".tr", self.nbins)):
            sums[cnt] = vals[:, :-1]
            nsnps_blk[cnt] = vals[:, -1].transpose()

        if not idx:
            self.nsnps_blk = nsnps_blk
        elif not (np.array_equal(self.nsnps_blk, nsnps_blk)):
            self.log._log("!!! Trace summary "+filename+" has different annotations !!!")
            sys.exit(1)
        
        self.sums.append(sums)
        self.ntrace += 1

        self.log._log(f"Read in trace summaries from {filename} generated with {K} random vectors")
        if (self.verbose):
            self.log._log("-- avg. jackknife LDscore sum:\n"+np.array2string(sums[:-1].mean(axis=0), precision=3, separator=', ')\
                +"\n-- number of jackknife blocks:\t"+str(self.nblks)\
                +"\n-- genome-wide LDscore sum:\n"+np.array2string(sums[-1], precision=3, separator=', '))
        return self.sums
        
    def _read_all_trace(self):
        trace_files = [self.sumpath]
        if path.isdir(self.sumpath):
            prefix = self.sumpath.rstrip('/') + "/"
            trace_files = sorted([prefix + f.rstrip('.tr') for f in listdir(self.sumpath) if f.endswith('.tr')])
        for i, f in enumerate(trace_files):
            self._read_trace(f, i)
        self.log._log("Finished reading "+str(len(trace_files))+" trace summaries.")
        if (np.std(self.sums, axis=0).sum() == .0) and (self.ntrace > 1):
            self.log._log("!!! Duplicate trace summaries are used -- effective number of random vectors remains unchanged. !!!")
            self.effective_K = np.min(self.K)
        else:
            self.effective_K = np.sum(self.K)
        ## TODO: is it correct to take a weighted average by sample size?
        self.sums = np.average(self.sums, axis=0, weights=self.nsamp)
        self.nsamp = np.mean(self.nsamp)
        return self.sums

    def _calc_trace(self, nsample):
        self.log._log("Calculating trace...")
        ## ensure that annotation dimension matches that of the trace summaries or LD scores
        if (self.ldscores is not None):
            return self._calc_trace_from_ldscores(nsample)
        else:
            return self._calc_trace_from_sums(nsample)
    
    def _calc_trace_from_sums(self, N):
        trace = np.full((self.nblks+1, self.nbins+1, self.nbins+1), N)
        for k in range(self.nbins):
            for l in range(self.nbins):
                for j in range(self.nblks+1):
                    trace[j, k, l] = utils._calc_trace_from_ld(self.sums[j, k, l], N, self.nsnps_blk[j, k], self.nsnps_blk[j, l])
        return trace

    def _reset(self):
        '''
        this is for running multiple phenotypes; reset the *_filt params
        '''
        self.nsnps_blk_filt = self.nsnps_blk
        self.sums_filt = self.sums
        self.nsnps_blk_filt = self.nsnps_blk

    def _filter_snps(self, removelist):
        '''
        Remove the SNPs in the removelist from trace calculation. If (truncated) LD scores are available,
        then use a scaled trace contribution from those SNPs; otherwise, assume uniform LD mapping (i.e.,
        long-distance LD distribution is the same for all the SNPs)
        Since different SNPs are filtered for each phenotype, modify only the "_filt" parameters
        '''
        # filter w/o ldscores
        if (self.ldscores is None):
            # filter w/o snplist (.bim): since we don't know which jn block
            # they belong to, treat as uniform chance of belonging to any block
            if (len(self.snplist) == 0):
                self.sums_filt *= (1 - len(removelist)/self.nsnps_blk_filt)
                self.nsnps_blk_filt -= len(removelist)*(self.nblks-1)/self.nblks
            # we know where the snp belongs to, but don't have their (truncated) LD information
            # then simply scale each block 
            else:
                self._calc_blk_ld()
                self._map_idx()
                removeidx = [self.mapping.get(snp, None) for snp in removelist]
                indices, counts = np.unique(removeidx, return_counts=True)
                # for idx, cnt in zip(indices, counts):
                #     self.ldsums_blk_filt[idx] *= (1 - cnt/self.nsnps_blk_filt[idx])
                #     self.nsnps_blk_filt[idx] -= cnt
                # print(sum(self.ldsums_blk_filt))
                # print(sum(self.nsnps_blk_filt))
                #updated_sums = self._calc_jn_subsample(self.ldsums_blk_filt)
                #print(updated_sums[:-1].min(), updated_sums[:-1].mean(), updated_sums.max(), updated_sums[-1])
        ### TODO: it's possible that some of the SNPs in the trace is not included in the sumstat
        ### in this case perhaps it's possible to remove those SNPs in the trace calculations as well
        ### implement this

    def _calc_trace_from_ldscores(self, N):
        trace = np.full((self.nblks+1, self.nbins+1, self.nbins+1), N)
        for k in range(self.nbins):
            for l in range(self.nbins):
                ld_sum = self.ldscores[self.annot[:, k]==1][:, l].sum()
                for j in range(self.nblks+1):
                    idx_start = self.blk_size*j
                    idx_end = self.nsnps if j==self.nblks-1 else self.blk_size*(j+1)
                    annot_jn = self.annot[idx_start:idx_end]
                    ldscores_jn = self.ldscores[idx_start:idx_end]
                    ld_sum_jn = ld_sum - ldscores_jn[annot_jn[:, k]==1][:, l].sum()
                    if (j == self.nblks):
                        ld_sum_jn = ld_sum
                    trace[j, k, l] = utils._calc_trace_from_ld(ld_sum_jn, N, self.nsnps_blk[j, k], self.nsnps_blk[j, l])
        return trace

    def _read_ldscores(self):
        '''
        Read the LD score matrix (X^T Xz) instead of trace summaries. Works with either the (truncated) LDSC LD scores (.l2.ldscore.gz) or
        the genome-wide LD scores (.gw.ldscore.gz)
        '''
        #self.ldscores_df = pd.read_csv(self.ldscorespath, compression='gzip', sep='\t', index_col=False)
        self.ldscores_df = pd.read_csv(self.ldscorespath, compression='gzip', delim_whitespace=True, index_col=False)
        self.ldscores = self.ldscores_df.iloc[:, 3:].to_numpy()
        self.snplist = self.ldscores_df['SNP'].to_numpy()
        self.nsnps = self.ldscores.shape[0]
        self.nbins = self.ldscores.shape[1]
        self.log._log("Loaded the LD score matrix with "+str(self.nsnps)+" SNPs and "+\
                        str(self.nbins)+" bins")
        

    # def _read_one_rhe(self, filename, idx):
    #     # read in metadata
    #     nsamp, nsnps, nblks, nbins = 0, 0, 0, 0
    #     # the number of blocks in all the trace stats must be the same
    #     with open(filename+".MN", 'r') as fd:
    #         next(fd)
    #         nsamp, nsnps, nblks, nbins = map(int, fd.readline().split(','))
    #     self.nsamp.append(nsamp)
    #     if not idx:
    #         self.nblks = nblks
    #         self.nbins = nbins
    #     elif (nblks != self.nblks):
    #         self.log._log("!!! RHE output "+filename+" has incorrect number of jackknife blocks !!!")
    #         return
    #     elif (nbins != self.nbins):
    #         self.log._log("!!! RHE output "+filename+" has incorrect number of annotation bins !!!")
    #         return
    #     trace = np.zeros((self.nblks+1, self.nbins, self.nbins))
    #     nsnps_blk = np.zeros((self.nblks+1, self.nbins))
        
    #     # read in trace values
    #     for cnt, vals in enumerate(utils._read_multiple_lines(filename+".trace", self.nbins)):
    #         trace[cnt] = vals[:, :-1]
    #         nsnps_blk[cnt] = vals[:, -1].transpose()
        
    #     # calculate sum of block LD for each jackknife subsample
    #     lsum = np.zeros((nblks+1, self.nbins, self.nbins))
    #     for k in range(self.nbins):
    #         for l in range(self.nbins):
    #             for j in range(self.nblks+1):
    #                 lsum[j, k, l] =  utils._calc_lsum(trace[j, k, l], nsamp, nsnps_blk[j, k], nsnps_blk[j, l])
    #     self.lsums.append(lsum)
    #     self.nrhe += 1
    #     # if idx is 0, save # of SNPs in each jackknife subsample (this should be the same in all outputs)
    #     if not idx:
    #         self.nsnps_blk = nsnps_blk
    #         self.nsnps = nsnps