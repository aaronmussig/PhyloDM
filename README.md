# 🌲 PhyloDM

[![PyPI](https://img.shields.io/pypi/v/phylodm?color=yellow)](https://pypi.org/project/phylodm/)
[![BioConda](https://img.shields.io/conda/vn/bioconda/phylodm?color=43b02a)](https://anaconda.org/bioconda/phylodm)
[![Crates](https://img.shields.io/crates/v/phylodm?color=orange)](https://crates.io/crates/phylodm)
[![DOI](https://zenodo.org/badge/251473194.svg)](https://zenodo.org/badge/latestdoi/251473194)

PhyloDM is a high-performance library that converts a phylogenetic tree into a pairwise distance matrix. 

For a tree with 30,000 taxa, PhyloDM will use:

* ~14GB of memory (94% less than DendroPy)
* ~1 minute of CPU time (183x faster than DendroPy).

PhyloDM is written in Rust and is exposed to Python via the Python PyO3 API. This means it 
can be used in either Python or Rust, however, the documentation below is written for use in Python. For Rust documentation, see [Crates.io](https://docs.rs/phylodm/latest/phylodm/).

## ⚙ Installation

*Requires Python 3.7+*

### PyPI

Pre-compiled binaries are packaged for most 64-bit Unix platforms. If you are installing on a different platform then you
will need to have [Rust](https://www.rust-lang.org/tools/install) installed to compile the binaries. 

```shell
python -m pip install phylodm
```

### Conda

```shell
conda install -c b bioconda phylodm
```


## 🐍 Quick-start

A pairwise distance matrix can be created from either a Newick file, or DendroPy tree.

```python
from phylodm import PhyloDM

# PREPARATION: Create a test tree
with open('/tmp/newick.tree', 'w') as fh:
    fh.write('(A:4,(B:3,C:4):1);')

# 1a. From a Newick file
pdm = PhyloDM.load_from_newick_path('/tmp/newick.tree')

# 1b. From a DendroPy tree
import dendropy
tree = dendropy.Tree.get_from_path('/tmp/newick.tree', schema='newick')
pdm = PhyloDM.load_from_dendropy(tree)

# 2. Calculate the PDM
dm = pdm.dm(norm=False)
labels = pdm.taxa()

"""
/------------[4]------------ A
+
|          /---------[3]--------- B
\---[1]---+
           \------------[4]------------- C
           
labels = ('A', 'B', 'C')
    dm = [[0. 8. 9.]
          [8. 0. 7.]
          [9. 7. 0.]]
"""
```

### Accessing data
The `dm` method generates a symmetrical NumPy matrix and returns a tuple of
keys in the matrix row/column order.

```python
# Calculate the PDM
dm = pdm.dm(norm=False)
labels = pdm.taxa()

"""
/------------[4]------------ A
+
|          /---------[3]--------- B
\---[1]---+
           \------------[4]------------- C
           
labels = ('A', 'B', 'C')
    dm = [[0. 8. 9.]
          [8. 0. 7.]
          [9. 7. 0.]]
"""

# e.g. The following commands (equivalent) get the distance between A and B
dm[0, 1]  # 8
dm[labels.index('A'), labels.index('B')]  # 8
```

### Normalisation

If the `norm` argument of `dm` is set to `True`, then the data will be normalised 
by the sum of all edges in the tree.


## ⏱ Performance
Tests were executed using `scripts/performance/Snakefile` on an Intel(R) Xeon(R) CPU E5-2650 v3 @ 2.30GHz.

For large numbers of taxa it is beneficial to use PhyloDM, however, if you have a small number 
of taxa in the tree it is beneficial to use DendroPy for the great features it provides.



![PhyloDM vs DendroPy resource usage](https://raw.githubusercontent.com/aaronmussig/PhyloDM/main/docs/img/performance.svg)
