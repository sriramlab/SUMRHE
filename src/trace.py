'''
Read the trace output from RHE & summarize into trace statistics

the RHE output metadata (.MN) should have the following:
    # individuals, # SNPs
    N, M
the RHE output tracedata (.trace) should have the following:
    jackknife subsample trace, # SNPs
    tr, M'
the SNP sets should be the same for all the outputs. 
it is recommended to input .bim for the list of SNPs (read in once)

the trace output (.sum
'''

import numpy as np
from os import listdir
from logger import Logger

def _calc_lsum(tr, n, m):
    return (tr - n)*pow(m,2)/pow(n,2)

class Trace:
    def __init__(self, bimpath = None, rhepath=None, sumpath=None, savepath=None, log=None):
        self.log = log
        self.rhepath=rhepath
        self.sumpath=sumpath
        self.savepath=savepath
        self.lsums = []
        self.nblks = 100
        self.nrhe = 0
        self.snplist = []
        self.nsnps_jn = []
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
        lsum = np.zeros(nblks)
        for i in range(nblks):
            lsum[i] =  _calc_lsum(trace[i], nsamp, nsnps_jn[i])
        
        self.lsums.append(lsum)
        self.nrhe += 1
        # if idx is 0, save # of SNPs in each jackknife subsample (this should be the same in all outputs)
        self.nsnps_jn = np.array(nsnps_jn)
        
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
            for i in range(self.nblks):
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
        self.nsnps_jn = np.array(nsnps_jn)
        self.nblks = len(self.sums)
        self.log._log("Reading in trace summaries from "+self.sumpath)
        self.log._log("-- mean lsum: "+format(self.sums.mean(), '.3f')+"\n-- N jackknife blocks: "+str(self.nblks))
        return self.sums, self.stds
    
    def _read_ldscore(self, path=None):
        ''' 
        TODO: for filtering SNP's. Adjust trace estimates based on (truncated) LD scores; this is optional
        '''
        if (path is None) or (path == ""):
            self.log._log("!!! Invalid LD score path given !!!")

    def _calc_trace(self, nsample):
        return self.sums * pow(nsample, 2) / pow(self.nsnps_jn, 2) + nsample
