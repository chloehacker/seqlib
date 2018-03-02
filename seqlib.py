import numpy as np
import pandas as pd

class seqlib:
    def __init__(self, ninds, nsites):
        self.ninds = ninds
        self.nsites = nsites
        self.seqs = self.simulate()
    
    def mutate(self, base):
        diff = set("ACTG") - set(base)
        return np.random.choice(list(diff))
        
    def simulate(self):
        ninds = self.ninds
        nsites = self.nsites
        oseq = np.random.choice(list("ACGT"), size=nsites) 
        arr = np.array([oseq for i in range(ninds)])   
        muts = np.random.binomial(1, 0.1, (ninds, nsites)) 
        for col in range(nsites):  
            newbase = mutate(arr[0, col]) 
            mask = muts[:, col].astype(bool) 
            arr[:, col][mask] = newbase 
        missing = np.random.binomial(1, 0.1, (ninds, nsites)) 
        arr[missing.astype(bool)] = "N"  
        return arr
    
    def filter_missing(self, maxfreq):   
        arr = self.seqs
        freqmissing = np.sum(arr == "N", axis=0) / arr.shape[0]  
        return arr[:, freqmissing <= maxfreq]
    
    def filter(self, minfreq, maxfreq):  
        maf = self.filter_maf(self.filter_missing(maxfreq), minfreq)
        return maf
    
    def filter_maf(self, minmaf):
        freqs = np.sum(arr != arr[0], axis=0) / arr.shape[0]
        maf = freqs.copy()
        maf[maf > 0.5] = 1 - maf[maf > 0.5]
        return arr[:, maf > minmaf]
    
    def maf(self):
        freqs = np.sum(arr != arr[0], axis=0) / arr.shape[0]
        maf = freqs.copy()
        maf[maf > 0.5] = 1 - maf[maf > 0.5]
        return freqs
    
    
    def calculcate_statistics(self):
        nd = np.var(arr == arr[0], axis=0).mean() #calculating the mean of the array
        mf = np.mean(np.sum(arr != arr[0], axis=0) / arr.shape[0]) #freq calculated by dividing sum of maf by whole sum
        inv = np.any(arr != arr[0], axis=0).sum() 
        var = arr.shape[1] - inv 
        return pd.Series(
            {"mean nucleotide diversity": nd,
             "mean minor allele frequency": mf,
             "invariant sites": inv,
             "variable sites": var,
            })      
