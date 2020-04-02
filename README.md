# PhyloDm

[![codecov.io](https://codecov.io/github/aaronmussig/PhyloDM/coverage.svg?branch=master)](https://codecov.io/github/aaronmussig/PhyloDM?branch=master)

Efficient calculation of pairwise phylogenetic distance matrices.

## Installation
* PyPI: `pip install phylodm`
* Conda: `conda install -c bioconda phylodm`

## Usage
The leaf nodes in the tree must have unique names, otherwise a `DuplicateIndex` exception is raised.

### Python library

#### Creating a PDM a DendroPy object
```python
import dendropy
from phylodm.pdm import PDM

tree = dendropy.Tree.get_from_string('(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);', 'newick')
```

`PDM.get_from_dendropy`


`PDM.get_from_newick_file`

`PDM.get_from_path`


### CLI
The CLI can be used to create a phylogenetic distance matrix given a newick tree, e.g.:
 
`python -m phylodm /path/to/newick.tree pd /path/to/matrix.mat`
