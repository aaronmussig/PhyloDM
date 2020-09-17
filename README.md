# PhyloDM
[![PyPI](https://img.shields.io/pypi/v/phylodm)](https://pypi.org/project/phylodm/)
![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/phylodm?color=green)
[![codecov.io](https://codecov.io/github/aaronmussig/PhyloDM/coverage.svg?branch=master)](https://codecov.io/github/aaronmussig/PhyloDM?branch=master)
[![DOI](https://zenodo.org/badge/251473194.svg)](https://zenodo.org/badge/latestdoi/251473194)

Efficient calculation of pairwise phylogenetic distance matrices.

## Installation
PhyloDM should work out of the box on *nix systems.

* PyPI: `pip install phylodm`
* conda: `conda install -c bioconda phylodm`

_Note: You must have a C compiler._

## Usage
The leaf nodes in the tree must have unique names, otherwise a `DuplicateIndex` exception is raised.

### Python library

#### Creating a phylogenetic distance matrix
A phylogenetic distance matrix (`PDM`) object can be created from either a DendroPy tree, or a
newick file:

```python
import dendropy
from phylodm.pdm import PDM

# Load from DendroPy
t = dendropy.Tree.get_from_string('(A:4,(B:3,C:4):1);', 'newick')
pdm = PDM.get_from_dendropy(tree=t, method='pd', cpus=1)

# Load from Newick
with open('/tmp/newick.tree', 'w') as fh:
    fh.write('(A:4,(B:3,C:4):1);')
pdm = PDM.get_from_newick_file('/tmp/newick.tree', method='pd', cpus=1)
```

Once created, a `PDM` can be cached to disk, where it can be later loaded:

```python
import dendropy
from phylodm.pdm import PDM

# Create a PDM.
t = dendropy.Tree.get_from_string('(A:4,(B:3,C:4):1);', 'newick')
pdm_a = PDM.get_from_dendropy(tree=t, method='pd', cpus=1)

# Write to cache.
pdm_a.save_to_path('/tmp/pdm.mat')

# Load from cache.
pdm_b = PDM.get_from_path('/tmp/pdm.mat')
```

#### Accessing data
The `PDM.as_matrix` method generates a symmetrical numpy distance matrix and returns a tuple of
keys in the matrix row/column order:
```python
import dendropy
from phylodm.pdm import PDM

# Load from DendroPy
t = dendropy.Tree.get_from_string('(A:4,(B:3,C:4):1);', 'newick')
pdm = PDM.get_from_dendropy(tree=t, method='pd', cpus=1)
labels, mat = pdm.as_matrix(normalised=False)
"""
/------------[4]------------ A
+
|          /---------[3]--------- B
\---[1]---+
           \------------[4]------------- C
           
labels = ('A', 'B', 'C')
mat = [[0. 8. 9.]
       [8. 0. 7.]
       [9. 7. 0.]]
"""

# Retrieving a specific value
pdm.get_value('A', 'A')  # 0
pdm.get_value('A', 'C')  # 9
pdm.get_value('C', 'A')  # 9
```

#### Method
The method parameter can be either patristic distance (`pd`) or the count of edges between 
leaves (`node`).

#### Normalisation
If true, the data will be returned as normalised depending on the method:
* `pd` = sum of all edges
* `node` = count of all edges

### CLI
The CLI can be used to create a phylogenetic distance matrix given a newick tree, e.g.:
 
`python -m phylodm /path/to/newick.tree pd /path/to/matrix.mat`

## Performance
Tests were executed using the `scripts/phylodm_perf.py` script with 10 trials.

These tests demonstrate that PhyloDM is more efficient than DendroPy's
phylogenetic distance matrix when there are over 500 taxa in the tree. If there
are less than 500 taxa, then use DendroPy for all of the great 
features it provides. 

With 10,000 taxa in the tree, each program uses approximately:
* PhyloDM = 4 seconds / 2 GB memory
* DendroPy = 17 minutes / 90 GB memory

![DendroPy vs. PhyloDM PDM Construction Time](docs/img/dendropy_vs_phylodm_time.png)![DendroPy vs. PhyloDM PDM Maximum Memory Usage](docs/img/denropy_vs_phylodm_memory.png)

## Changelog
```
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
