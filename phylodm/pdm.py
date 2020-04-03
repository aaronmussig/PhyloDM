###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

from collections import deque
from typing import Tuple, Union

import dendropy
import h5py
import numpy as np
from tqdm import tqdm

from phylodm.common import create_mat_vector, compact_int_mat
from phylodm.indices import Indices
from phylodm.symmat import SymMat


class PDM(SymMat):

    def __init__(self):
        """Variables are set when called by any of the static methods."""
        super().__init__()
        self._tree_length = None
        self._method = None

    @staticmethod
    def get_from_dendropy(tree: dendropy.Tree, method: str = 'pd') -> 'PDM':
        """Create a PDM from a dendropy tree object."""
        return PDM()._get_from_dendropy(tree, method)

    @staticmethod
    def get_from_newick_file(path_tree: str, method: str = 'pd') -> 'PDM':
        """Create a PDM from a newick file."""
        tree = dendropy.Tree.get_from_path(path_tree,
                                           schema='newick',
                                           rooting='force-unrooted',
                                           preserve_underscores=True)
        return PDM.get_from_dendropy(tree, method)

    @staticmethod
    def get_from_path(path: str) -> 'PDM':
        """Load a PDM from a cache."""
        return PDM()._get_from_path(path)

    def _get_from_dendropy(self, tree: dendropy.Tree, method: str = 'pd') -> 'PDM':
        self._method = method
        if method == 'pd':
            use_pd = True
            self._d_type = np.dtype('float64')
            self._arr_default = 0.0
        elif method == 'node':
            use_pd = False
            self._d_type = np.dtype('uint32')
            self._arr_default = 0
        else:
            raise NotImplemented(f'Unknown method: {method}')

        # Deep copy tree
        tree = dendropy.Tree.get_from_string(tree.as_string(schema='newick'), 'newick')

        # Check that the end condition can be satisfied.
        if tree.seed_node.parent_node is not None:
            raise Exception('Dendropy did not seed the tree correctly.')

        # Store the tree length for normalising branch lengths.
        if use_pd:
            self._tree_length = tree.length()
        else:
            self._tree_length = len(tree.edges())

        # Queue each of the leaf nodes for processing.
        queue = deque()
        for leaf_node in tree.leaf_node_iter():
            queue.append(leaf_node)

        # Create leaf indices
        self._indices = Indices()
        for leaf in sorted(queue, key=lambda x: x.taxon.label):
            self._indices.add_key(leaf.taxon.label)

        self._data = create_mat_vector(len(self._indices), self._arr_default)

        # Track distances until the seed node, from each leaf node.
        for leaf_node in tqdm(queue):
            cur_node = leaf_node
            cur_node.child_dist = {leaf_node: 0.0 if use_pd else 0}
            sisters_processed = set()

            while cur_node.parent_node is not None:
                parent_node = cur_node.parent_node
                tot_dist = cur_node.child_dist[leaf_node] + (cur_node.edge_length if use_pd else 1)

                # First visit, just add the cumulative distance.
                if not hasattr(parent_node, 'child_dist'):
                    parent_node.child_dist = {}

                # Calculate all pairwise distances between sister nodes not processed
                for sister_node in set(parent_node.child_dist.keys()).difference(sisters_processed):
                    sister_dist = parent_node.child_dist[sister_node]

                    if leaf_node.taxon.label == sister_node.taxon.label:
                        calc_dist = self._arr_default
                    else:
                        calc_dist = sister_dist + tot_dist
                    self.set_value(leaf_node.taxon.label, sister_node.taxon.label, calc_dist)

                    sisters_processed.add(sister_node)

                parent_node.child_dist[leaf_node] = tot_dist
                cur_node = parent_node

        # Use the smallest possible data type for the matrix
        if method == 'node':
            self._data = compact_int_mat(self._data)
            self._d_type = self._data.dtype
        return self

    def as_matrix(self, normalised: bool = False) -> Tuple[Tuple[str], np.array]:
        """Return the PDM as a symmetric numpy matrix."""
        labels, mat = super().as_matrix()
        if normalised:
            mat = mat * (1.0 / self._tree_length)
        return labels, mat

    def save_to_path(self, path: str):
        """Write the PDM to a cache."""
        with h5py.File(path, 'w') as f:
            f.create_dataset('indices',
                             data=[t.encode('ascii') for t in self._indices.get_keys()],
                             dtype=h5py.string_dtype(encoding='ascii'))
            f.create_dataset('data', data=self._data, dtype=self._d_type)
            f.create_dataset('tree_length', data=self._tree_length)
            f.create_dataset('arr_default', data=self._arr_default)
            f.create_dataset('method', data=self._method.encode('ascii'), dtype=h5py.string_dtype(encoding='ascii'))

    def _get_from_path(self, path: str) -> 'PDM':
        """Read the PDM from a cache."""
        self._indices = Indices()
        with h5py.File(path, 'r') as hf:
            self._arr_default = hf['arr_default'][()]
            self._data = hf['data'][()]
            self._tree_length = hf['tree_length'][()]
            self._method = hf['method'][()].decode('utf-8')
            for idx, key in enumerate(hf['indices'][()]):
                self._indices.add_key(key.decode('utf-8'))
        self._d_type = self._data.dtype
        return self

    def get_value(self, key_i: str, key_j: str, normalised: bool = False) -> Union[int, float]:
        """Get the value of two keys."""
        value = super().get_value(key_i, key_j)
        if normalised:
            value *= 1 / self._tree_length
        return value
