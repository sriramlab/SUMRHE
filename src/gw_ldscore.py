"""
Stochastically estimate (partitioned) genome-wide LD scores. Part of the code is modified from Eric Liu's script
"""
import utils

import numpy as np
import pandas as pd
from bed_reader import open_bed
import multiprocessing as mp
from tqdm import tqdm
import scipy
import sys
import gc

class GenomewideLDScore:
    def __init__(self, bed_path, annot_path, out_path, log, num_vecs=10, num_workers=4, step_size=1000):
        self.G = open_bed(bed_path+".bed")
        self.nsamp, self.nsnps = self.G.shape
        self.nvecs = num_vecs
        self.nworkers = num_workers
        self.step_size = step_size
        self.snp_idx = np.arange(self.nsnps)
        self.log = log
        self._read_bim(bed_path+".bim")
        self.log._log("Calculating genome-wide LD scores from genotype file: "+str(bed_path))
        self.log._log("Number of individuals: "+str(self.nsamp)+", Number of SNPs: "+str(self.nsnps))
        self.outpath = out_path
        if (annot_path is not None):
            self.log._log("With annotation file: "+str(annot_path))
        self._read_annot(annot_path)

    def _compute_Xz_blk(self, blk_idxs):
        """
        Compute Xk z given partitioned indices for the blk. Returns (N x nblks x nvecs)
        Read in the genotype blk altogether, then splice it afterwards (to make things faster).
        idx passed is the "blk" idx, i.e., always starts from 0.
        """
        blk_start, blk_end, idxs = blk_idxs # j == blk idx, idxs == nested list of snp idx
        nsnps = sum([len(binidx) for binidx in idxs])
        Xz = np.zeros((self.nbins, self.nsamp, self.nvecs))
        Zs = np.random.normal(size=(nsnps, self.nvecs))

        geno = self.G.read(index=np.s_[:, blk_start: blk_end])
        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)

        geno = (geno-means)/stds
        geno[np.isnan(geno)] = 0

        geno = np.array(geno, order='F')
        Zs = np.array(Zs, order='F')
        
        for k, binidx in enumerate(idxs):
            Xz[k, :, :] = scipy.linalg.blas.sgemm(1.0, geno[:, binidx], Zs[binidx, :])
        return Xz
    

    def _compute_XtXz_blk(self, blk_idx):
        """
        For blk genotype, multiply with X_k z to get XtXkz.
        """
        blk_start, blk_end = blk_idx

        geno = self.G.read(index=np.s_[:, blk_start:blk_end])
        means = np.nanmean(geno, axis=0)
        stds = np.nanstd(geno, axis=0)

        geno = (geno-means)/stds
        geno[np.isnan(geno)] = 0

        ## TODO: benchmark sgemm vs. np broadcasting - in small scale, looks like sgemm is faster (could be b/c sgemm is used in Xz estimation)
        geno_t = np.array(geno.T, order='F')
        XtXz = np.zeros((blk_end - blk_start, self.nbins, self.nvecs))
        for k in range(self.nbins):
            XtXz[:, k, :] = scipy.linalg.blas.sgemm(1.0, geno_t, self.Xz[k])

        #XtXz = np.einsum('nm,knb->mkb', geno, self.Xz)
        return (blk_start, blk_end, XtXz)


    def _read_annot(self, annot_path):
        """
        Read in the annotation. If the file includes a header, save it as the names for the annotations.
        If not, then have dummy names and read in the annotation.
        """
        if (annot_path is None):
            self.annot = np.ones((self.nsnps, 1))
            self.log._log("Calculating genome-wide (non-partitioned) LD score")
        else:
            self.l2cols, self.annot = utils._read_with_optional_header(annot_path)
            if (self.annot.ndim == 1):
                self.annot = self.annot.reshape(-1, 1)
            self.log._log("Read SNP partition annotation of dimensions "+str(self.annot.shape))
        
        if (self.nsnps != self.annot.shape[0]):
            self.log._log(f"!!! number of SNPs in annotation ({self.annot.shape[0]}) does not match the input genotype file ({self.nsnps}) !!!")
            sys.exit(1)
        self.nbins = self.annot.shape[1]
        if (self.l2cols is None):
            self.l2cols = ['L2_'+str(i) for i in range(self.annot.shape[1])]
        else:
            self.l2cols = [i + 'L2' for i in self.l2cols]
        self.log._log(f"Nbins: {self.nbins}")
        self.nsnps_bin = self.annot.sum(axis=0)

    def _read_bim(self, bim_path):
        if (bim_path is None):
            self.log._log("No .bim file is provided; all (anonymous) SNPs will be used")
            self.snplist = None
        else:
            self.log._log(f"Reading {bim_path} for SNPs")
            self.snplist = pd.read_csv(bim_path, header=None, sep='\t')
        if (len(self.snplist) != self.nsnps):
            self.log._log(f"!!! The number of SNPs in the .bed file ({self.nsnps}) does not match the .bim file ({len(self.snplist)}) !!!")
            sys.exit(1)
    
    def _partition_index(self, snpidx, annot) -> list[np.ndarray]:
        """
        partition snp indices by annotation
        """
        return [snpidx[annot[:, c] == 1] for c in range(self.nbins)]


    def _compute_ldscore(self):
        """
        Use multi-processing to calculate the X_j^T X_k Z.
        General sketch: read in each block of genotype, calculate X_k Z for that blk. Aggregate X_k through all blks.
        Then re-read each blk from the start, multiply by the previous result (loop over k) to get X_j ^ T X_k (no need for agg this time).
        """
        self.start_time = utils._get_time()
        self.log._log("Genome-wide LD score calculation started at: "+utils._get_timestr(self.start_time))
        self.log._log(f"num_vecs: {self.nvecs}, num_workers: {self.nworkers}, step_size: {self.step_size}")

        self.nblks = len(np.arange(self.nsnps)[::self.step_size])
        Xz_input = []
        XtXz_input = []
        for j in range(self.nblks):
            idx_start = self.step_size*j
            idx_end = self.nsnps if j==self.nblks-1 else self.step_size*(j+1)
            annot_blk = self.annot[idx_start:idx_end]
            Xz_input.append((idx_start, idx_end, self._partition_index(np.arange(len(annot_blk)), annot_blk)))
            XtXz_input.append((idx_start, idx_end))
        
        self.Xz = np.zeros((self.nbins, self.nsamp, self.nvecs))

        with mp.Pool(self.nworkers) as pool:
            with tqdm(total=self.nblks) as pbar:
                pbar.set_description('Calculating Xz')
                for result in pool.imap_unordered(self._compute_Xz_blk, Xz_input):
                    self.Xz += result
                    gc.collect()
                    pbar.update()
        pool.join()

        self.Xz_time = utils._get_time()
        self.log._log("Calculation of Xz (for each partition) completed. Runtime: "+format(self.Xz_time - self.start_time, '.3f')+" s")

        self.XtXz = np.zeros((self.nsnps, self.nbins, self.nvecs))

        with mp.Pool(self.nworkers) as pool:
            with tqdm(total=self.nblks) as pbar:
                pbar.set_description('Calculating XtXz')
                for result in pool.imap_unordered(self._compute_XtXz_blk, XtXz_input):
                    idx_start, idx_end, XtXz_blk = result
                    self.XtXz[idx_start:idx_end, :, :] = XtXz_blk
                    gc.collect()
                    pbar.update()
        pool.join()

        self.XtXz = self.XtXz / self.nsamp

        self.XtXz_time = utils._get_time()
        self.log._log("Calculation of XtXz (for each partition) completed. Runtime: "+format(self.XtXz_time - self.Xz_time, '.3f')+" s")

        self.log._log("Converting XtXz into genome-wide (partitioned) LD scores.")
        self.gwldscore = self.nsamp/(self.nsamp+1) * (np.square(self.XtXz).mean(axis=2) - self.nsnps_bin/self.nsamp)
        self.log._log(f"Saving the genome-wide (partitioned) LD scores into: {self.outpath}.gw.ldscore.gz")
        snpcols = ['CHR', 'SNP', 'BP']
        if (self.snplist is None):
            self.snpdf = pd.DataFrame(np.nan*np.ones((self.nsnps, 3)), columns=snpcols)
        else:
            self.snpdf = self.snplist.iloc[:, :3]
            self.snpdf.columns = snpcols
        
        self.gwldscore = pd.DataFrame(self.gwldscore, columns = self.l2cols)
        self.gwldscore = pd.concat([self.snpdf, self.gwldscore], axis=1)
        self.gwldscore.to_csv(f'{self.outpath}.gw.ldscore.gz', index=False, compression='gzip', sep='\t', float_format='%.3f')
        self.end_time = utils._get_time()
        self.log._log(f"Calculation of genome-wide (non-partitioned) LD score ended at "+utils._get_timestr(self.end_time))
        self.log._log("Runtime: "+format(self.end_time - self.start_time, '.3f')+" s")
        self.log._save_log(self.outpath+".log")

        





        



        

            