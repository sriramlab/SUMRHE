import numpy as np

def _partition_bin_non_overlapping(jn_values: np.ndarray, jn_annot: np.ndarray, nbins: int):
    """
    partition the first array (a 1D np array) by the annotation (a 1D array). return a nested list.
    This function assumes that a SNP can belong to only one bin.
    """
    partitions = {i: [] for i in range(nbins)}
    for z, idx in zip(jn_values, jn_annot):
        partitions[idx.argmax()].append(z)
    return [partitions[i] for i in range(nbins)], [len(partitions[i]) for i in range(nbins)]


def _partition_bin_overlapping(jn_values: np.ndarray, jn_annot: np.ndarray, nbins: int):
    """
    partition the first array (a 1D np array) by the annotation (a 1D array). return a nested list.
    This function assumes that a SNP can belong to multiple bins.
    """
    partitions = {i: [] for i in range(nbins)}
    for bin in range(nbins):
        partitions[bin] = [jn_values[snp] for snp in range(len(jn_values)) if jn_annot[row, bin]]
    return [partitions[i] for i in range(nbins)], [len(partitions[i]) for i in range(nbins)]


# def _partition_count(jn_annot: np.ndarray, nbins: int):
#     """
#     return the partitioned count
#     """
#     counts = np.zeros(nbins, dtype=int)
#     value_cnts = 
