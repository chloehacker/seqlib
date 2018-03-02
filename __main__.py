class Seqlib:
    def __init__(self, ninds, nsites, arr):
        self.ninds = ninds
        self.nsites = nsites
        self.seqs = self.simulate()
        self.arr = arr
        
    def simulate(self):
        pass
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
    
    def mutate(base):
        diff = set("ACTG") - set(base)
        return np.random.choice(list(diff))
   
    def filter_missing(self, arr, maxfreq):   
        freqmissing = np.sum(arr == "N", axis=0) / arr.shape[0]  
        return arr[:, freqmissing <= maxfreq]
    
    def filter_maf_missing(self,minfreq, maxfreq):     
        return filter_maf(filter_missing(maxfreq), minfreq)
    
    
    def calculcate_statistics(arr):
        nd = np.var(arr == arr[0], axis=0).mean() #calculating the mean of the array
        mf = np.mean(np.sum(arr != arr[0], axis=0) / arr.shape[0]) #freq calculated by dividing sum of maf by whole sum
        inv = np.any(arr != arr[0], axis=0).sum() # sums invariant sites
        var = arr.shape[1] - inv #subtracts invariant from total
        return pd.Series(
            {"mean nucleotide diversity": nd,
             "mean minor allele frequency": mf,
             "invariant sites": inv,
             "variable sites": var,
            })      