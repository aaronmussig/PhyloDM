# PhyloDM
[![PyPI](https://img.shields.io/pypi/v/phylodm)](https://pypi.org/project/phylodm/)
[![BioConda](https://img.shields.io/conda/vn/bioconda/phylodm?color=green)](https://anaconda.org/bioconda/phylodm)
[![Crates](https://img.shields.io/crates/v/phylodm?color=orange)](https://crates.io/crates/phylodm)
[![DOI](https://zenodo.org/badge/251473194.svg)](https://zenodo.org/badge/latestdoi/251473194)

Efficient calculation of pairwise phylogenetic distance matrices.

## Installation

The easiest installation method is through Conda. 
If you choose to install via PyPI 
ensure that you have a [Rust compiler](https://www.rust-lang.org/tools/install).

* PyPI: `pip install phylodm`
* Conda: `conda install -c bioconda phylodm`


## Usage

### Creating a phylogenetic distance matrix
A phylogenetic distance matrix (`PhyloDM`) object can be created from a newick file:

```python
from phylodm import PhyloDM

# Create a test tree
with open('/tmp/newick.tree', 'w') as fh:
    fh.write('(A:4,(B:3,C:4):1);')

# Load from newick
pdm = PhyloDM()
pdm.load_from_newick_path('/tmp/newick.tree')
```

### Accessing data
The `dm` method generates a symmetrical numpy distance matrix and returns a tuple of
keys in the matrix row/column order:

```python
from phylodm import PhyloDM

# Create a test tree
with open('/tmp/newick.tree', 'w') as fh:
    fh.write('(A:4,(B:3,C:4):1);')

# Load from Newick file
pdm = PhyloDM.load_from_newick_path('/tmp/newick.tree')

# Or, load from a Dendropy object
import dendropy
tree = dendropy.Tree.get_from_path('/tmp/newick.tree', schema='newick')
pdm = PhyloDM.load_from_dendropy(tree)

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
```


### Normalisation
If true, the data will be returned as normalised by the sum of all edges in the tree.


## Performance
Tests were executed using the `scripts/phylodm_perf.py` script with 10 trials.

These tests demonstrate that PhyloDM is more efficient than DendroPy's
phylogenetic distance matrix when there are over 500 taxa in the tree. If there
are less than 500 taxa, then use DendroPy for all of the great 
features it provides. 

With 10,000 taxa in the tree, each program uses approximately:
* PhyloDM = 4 seconds / 40 MB memory
* DendroPy = 17 minutes / 22 GB memory

![DendroPy vs. PhyloDM PDM Construction Time](docs/img/dendropy_vs_phylodm_time.png)![DendroPy vs. PhyloDM PDM Maximum Memory Usage](docs/img/denropy_vs_phylodm_memory.png)

## Changelog
```
2.0.0
  - Re-write in Rust (2x faster)
1.3.1
  - Use OpenMP to parallelize PDM methods.
1.3.0
  - Removed tqdm.
  - get_matrix() is now 3x faster.
1.2.0
  - Addded the remove_keys command.
1.1.0
  - Significant improvement in PDM construction time using C.
1.0.0
  - Initial release.
```

## Citing
Please [cite this software](https://doi.org/10.5281/zenodo.3998716) if you use it in your work.
