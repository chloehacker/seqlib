# seqlib library
seqlib class object assignment

### To install:

```
git clone https://github.com/chloehacker/seqlib.git
cd seqlib
pip install .
```

### To use:

```
# import the library
import seqlib

# class object generate with given shape
s = seqlib.seqlib(10, 100)

# prints sequence array
print(s.seqs)

# prints minor allele frequency of sequence array
print (s.maf)

# prints filtered array based on map and missing sites
print(s.filter(minmaf=0.1, maxmissing=0.0))

# new copy of filtered array
n = s.filter_seqlib(minmaf=0.1, maxmissing=0.0)

# stats of full array
s.calculate_statistics()

# stats of filtered array
n.calculate_statistics()

# both in one
s.filter_seqlib(minmaf=0.1, maxmissing=0.0).calculate_statistics()
```

