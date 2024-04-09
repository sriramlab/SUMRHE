'''
Read the trace output from RHE & summarize into trace statistics

the RHE output metadata (.MN) should have the following:
    # individuals, # SNPs
    N, M
the RHE output tracedata (.trace) should have the following:
    jackknife subsample trace, # SNPs
    tr, M'
There should be (njackknife+1) rows excluding the header, where
the last row is the trace with the entire genotype.

The SNP sets should be the same for all the outputs. 
it is recommended to input .bim for the list of SNPs (read in once)
'''

import numpy as np
from os import listdir
from logger import Logger

def _calc_lsum(tr, n, m):
    return (tr - n)*pow(m,2)/pow(n,2)

class Trace:
    def __init__(self, bimpath = None, rhepath=None, sumpath=None, savepath=None, log=None, ldproj=None, nblks=100):
        self.log = log
        self.rhepath = rhepath
        self.sumpath = sumpath
        self.savepath = savepath
        self.ldprojpath = ldproj
        self.lsums = []
        self.sums = None
        self.nblks = nblks # nblks specified only if using ld proj; otherwise it'll be overwritten by trace summaries.
        self.nrhe = 0
        self.snplist = []
        self.nsnps = 0
        self.nsnps_jn = None # this is the nsnps for each jn. does not change
        self.nsnps_jn_filt = None # this is the nsnps for each jn, which is filtered in case of --filter-both-sides
        self.ldscores = None
        if (bimpath is None) or (bimpath == ""):
            self.log._log("!!! SNP list (.bim) is generally recommended !!!")
        elif (bimpath.endswith(".bim")):
            with open(bimpath, 'r') as fd:
                for line in fd:
                    self.snplist.append(line.split()[1])
        else:
            self.log._log('!!! '+bimpath+' is not a .bim file! !!!')

    def _read_one_rhe(self, filename, idx):
        trace = []
        nsnps_jn = []
        nsamp = 0
        nsnps = 0
        nblks = 0 # at the moment, the number of blocks in all the trace stats must be the same
        # read in metadata
        with open(filename+".MN", 'r') as fd:
            next(fd)           
            nsamp, nsnps, nblks = map(int, fd.readline().split(','))
        if not idx:
            self.nblks = nblks
        elif (nblks != self.nblks):
            self.log._log("!!! RHE output "+filename+" has incorrect number of jackknife blocks !!!")
            return
        
        # read in trace values
        with open(filename+".trace", 'r') as fd:
            next(fd)
            for line in fd:
                tr, ns = map(float, line.strip().split(','))
                trace.append(tr)
                nsnps_jn.append(ns)
        
        # calculate sum of block LD for each jackknife subsample
        lsum = np.zeros(nblks+1)
        for i in range(nblks+1):
            lsum[i] =  _calc_lsum(trace[i], nsamp, nsnps_jn[i])
        
        self.lsums.append(lsum)
        self.nrhe += 1
        # if idx is 0, save # of SNPs in each jackknife subsample (this should be the same in all outputs)
        self.nsnps_jn = np.array(nsnps_jn)
        self.nsnps = nsnps
        
    def _read_all_rhe(self, path=None):
        if self.rhepath is None:
            self.rhepath = path
        trace_files = sorted([f.rstrip('.trace') for f in listdir(self.rhepath) if f.endswith('.trace')])
        prefix = self.rhepath.rstrip('/') + "/"
        for i, f in enumerate(trace_files):
            self._read_one_rhe(prefix+f, i)
        self.lsums = np.array(self.lsums)
        self.log._log("Finished reading "+str(len(trace_files))+" RHE trace outputs.")
        self.sums = self.lsums.mean(axis=0)
        self.stds = self.lsums.std(axis=0)
        return self.sums, self.stds
    
    def _save_trace(self):
        ''' Save trace summaries as a file'''
        with open(self.savepath+".tr", 'w') as fd:
            fd.write("LDsum,LDsum_std,NSNP\n")
            for i in range(self.nblks+1):
                fd.write(format(self.sums[i],'.3f')+","+format(self.stds[i], '.3f')+","+format(self.nsnps_jn[i], '.0f')+"\n")
        self.log._log("Saved trace summary into "+self.savepath+".tr")
    
    def _read_trace(self):
        ''' read from a trace summary file'''
        sums = []
        stds = []
        nsnps_jn = []
        with open(self.sumpath , 'r') as fd:
            next(fd)
            for line in fd:
                s, st, nj = map(float, line.strip().split(','))
                sums.append(s)
                stds.append(st)
                nsnps_jn.append(nj)
                
        self.sums = np.array(sums)
        self.stds = np.array(stds)
        self.nblks = len(self.sums)-1
        self.nsnps = nsnps_jn[self.nblks]
        self.nsnps_jn = np.array(nsnps_jn)
        self.log._log("Reading in trace summaries from "+self.sumpath)
        self.log._log("-- avg. jackknife LDscore sum:\t"+format(self.sums[:-1].mean(), '.3f')+"\n-- number of jackknife blocks:\t"+str(self.nblks)\
            +"\n-- genome-wide LDscore sum:\t"+format(self.sums[-1], '.3f'))
        return self.sums, self.stds
    
    def _read_ldscore(self, path=None):
        ''' 
        For filtering SNP's. Adjust trace estimates based on (truncated) LD scores; this is optional
        LD scores should be in the standard LD score format (CHR SNP BP L2)
        '''
        if (path == ""):
            self.log._log("!!! Invalid LD score path given. Proceeding filtering without ld scores !!!")
            return
        else:
            ldscores = []
            with open(path, 'r') as fd:
                fd.next()
                for line in fd:
                    ldscores.append(line.split()[3])
            # save ldscores as a dictionary
            self.ldscores = dict(zip(self.snplist, ldscores))
        return
            

    def _calc_trace(self, nsample):
        self.log._log("Calculating trace...")
        if (self.ld_proj is not None):
            return self._calc_trace_from_ldproj(nsample)
        else:
            #return self.sums * pow(nsample, 2) / pow(self.nsnps_jn_filt, 2) + nsample
            return self.sums_filt * pow(nsample, 2) / pow(self.nsnps_jn_filt, 2) + nsample

    @staticmethod
    def _partition(alist, npartition):
        step_size = len(alist)//npartition
        return [alist[i*step_size:(i+1)*step_size] if i+1 < npartition else alist[i*step_size: ] for i in range(npartition)]
        
    def _map_idx(self):
        '''
        create a mapping of SNP id -> idx
        '''
        self.mapping = {}
        partition = self._partition(self.snplist, self.nblks)
        for idx, part in enumerate(partition):
            for snp in part:
                self.mapping[snp] = idx
        return self.mapping

    @staticmethod
    def _calc_jn_subsample(alist):
        '''
        from a list/array return an array of leave-one-out (jackknife) subsamples
        the last element is the sum of all elements
        '''
        total = sum(alist)
        jn_sub = [total - val for val in alist]
        jn_sub.append(total)
        return np.array(jn_sub)

    # def _calc_blk_ld(self):
    #     total_ldscore = self.sums[self.nblks]
    #     self.ldsums_blk_filt = [total_ldscore - jack_ldscore for jack_ldscore in self.sums[:-1]]
    #     self.nsnps_blk_filt = np.array([self.nsnps_jn[self.nblks] - nsnp for nsnp in self.nsnps_jn[:-1]])
    #     return

    def _reset(self):
        '''
        this is for running multiple phenotypes; reset the *_filt params
        '''
        self.nsnps_jn_filt = self.nsnps_jn
        self.sums_filt = self.sums
        self.nsnps_jn_filt = self.nsnps_jn

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
                self.sums_filt *= (1 - len(removelist)/self.nsnps_jn_filt)
                self.nsnps_jn_filt -= len(removelist)*(self.nblks-1)/self.nblks
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

    @staticmethod
    def _calc_trace_from_ld(mat, b, n, m):
        #return np.square(mat).sum()*pow(n, 2)/(b*pow(m,2))
        return mat.sum()*pow(n/m, 2)/b + n

    def _calc_trace_from_ldproj(self, N):
        trace = np.zeros(self.nblks+1)
        for j in range(self.nblks):
            mask = np.ones((self.nsnps, ), dtype=bool)
            blk_size = self.nsnps//self.nblks
            idx_start = blk_size*j
            idx_end = 0
            if (j==self.nblks-1):
                idx_end = self.nsnps
            else:
                idx_end = blk_size*(j+1)
            mask[idx_start:idx_end] = False
            trace[j] = self._calc_trace_from_ld(self.ld_proj[mask], self.B, N, self.nsnps_jn_filt[j])
        trace[self.nblks] = self._calc_trace_from_ld(self.ld_proj, self.B, N, self.nsnps)
        return trace

    def _read_ldproj(self):
        '''
        Read the LD projection matrix (X^T Xz) instead of trace summaries. LD projection matrices
        are in np binary format (.npy)
        '''
        self.ld_proj = np.load(self.ldprojpath, allow_pickle=False)
        self.nsnps = self.ld_proj.shape[0]
        self.B = self.ld_proj.shape[1]
        self.log._log("Loaded the LD projection matrix with "+str(self.nsnps)+" SNPs and "+\
                        str(self.B)+" random vectors")
        blk_size = self.nsnps//self.nblks
        self.log._log("blk_size: "+str(blk_size))
        self.nsnps_jn = [self.nsnps - blk_size if i < self.nblks-1 else i*blk_size for i in range(self.nblks)]
        self.nsnps_jn.append(self.nsnps)
        return
