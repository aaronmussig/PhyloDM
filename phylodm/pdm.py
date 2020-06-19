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

from collections import deque, defaultdict
from typing import Tuple, Union

import dendropy
import h5py
import numpy as np
from tqdm import tqdm

from phylodm.common import create_mat_vector, compact_int_mat
from phylodm.indices import Indices
from phylodm.pdm_c import cartesian_sum
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

    @staticmethod
    def _get_int_node_depth(tree):
        out = defaultdict(set)
        queue = deque([(0, tree.seed_node)])
        while len(queue) > 0:
            depth, node = queue.pop()
            out[depth].add(node)
            for child_node in node.child_node_iter():
                if not child_node.is_leaf():
                    queue.append((depth + 1, child_node))
        return out

    @staticmethod
    def _preprocess_tree(tree: dendropy.Tree) -> dendropy.Tree:
        """Deep copy and pre-process a DendroPy tree."""
        out = dendropy.Tree.get_from_string(tree.as_string(schema='newick'), 'newick')
        if tree.seed_node.parent_node is not None:
            raise Exception('DendroPy did not seed the tree correctly.')
        return out

    def _get_from_dendropy(self, tree: dendropy.Tree, method: str = 'pd') -> 'PDM':

        # Initialise the tree
        tree = PDM._preprocess_tree(tree)

        # Initialise the class based on arguments.
        self._method = method
        if method == 'pd':
            use_pd = True
            self._d_type = np.dtype('float64')
            self._arr_default = 0.0
            self._tree_length = tree.length()
        elif method == 'node':
            use_pd = False
            self._d_type = np.dtype('uint32')
            self._arr_default = 0
            self._tree_length = len(tree.edges())
        else:
            raise NotImplemented(f'Unknown method: {method}')

        # Get each node at a specific depth.
        depth_to_node = PDM._get_int_node_depth(tree)

        # Pre-process each leaf node and create the indices.
        self._indices = Indices()
        for leaf_node in sorted(tree.leaf_node_iter(), key=lambda x: x.taxon.label):
            leaf_idx = self._indices.add_key(leaf_node.taxon.label)
            leaf_node.attr_child_dist = [(leaf_idx, leaf_node.edge_length if use_pd else 1)]

        # Create the vector to store the results.
        self._data = create_mat_vector(len(self._indices), 0.0)

        # Process the deepest nodes first, merging data at each level.
        with tqdm(total=len(self._data), unit_scale=True) as p_bar:
            n_indices = len(self._indices)
            for cur_depth, cur_nodes in sorted(depth_to_node.items(), key=lambda x: -x[0]):
                for node in cur_nodes:

                    # Extract the groups
                    groupings = [0]
                    group_idxs = list()
                    group_vals = list()
                    child_dist, child_groups = list(), list()
                    offset = 0

                    child_groups = list()

                    for child_node in node.child_node_iter():
                        child_groups.append(child_node.attr_child_dist)
                        for desc_idx, desc_dist in child_node.attr_child_dist:
                            group_idxs.append(desc_idx)
                            group_vals.append(desc_dist)
                            offset += 1
                        groupings.append(offset)
                        child_dist.extend(child_node.attr_child_dist)
                        del child_node.attr_child_dist

                    groupings = np.array(groupings, dtype=np.uint32)
                    group_idxs = np.array(group_idxs, dtype=np.uint32)
                    group_vals = np.array(group_vals, dtype=np.float64)

                    n_iter = cartesian_sum(n_indices, groupings, group_idxs, group_vals, self._data)
                    p_bar.update(n_iter)

                    # Record the distance from this node to its leaf nodes, and bring up.
                    # Add th distance from thsi node to its children
                    if cur_depth > 0:
                        node.attr_child_dist = [(x, y + (node.edge_length if use_pd else 1)) for x, y in child_dist]

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
