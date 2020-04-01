from collections import deque
from typing import Tuple, List

import dendropy
import h5py
import numpy as np
from tqdm import tqdm

from phylodm.common import create_mat
from phylodm.indices import Indices
from phylodm.symmat import SymMat


class PDM(SymMat):

    def __init__(self):
        super().__init__()
        self._tree_length = None

    def as_matrix(self, normalised: bool = False) -> Tuple[List[str], np.array]:
        labels, mat = super().as_matrix()
        if normalised:
            mat = mat * (1.0 / self._tree_length)
        return labels, mat

    @staticmethod
    def get_from_dendropy(tree: dendropy.Tree, method: str = 'pd') -> 'PDM':
        return PDM()._get_from_dendropy(tree, method)

    def _get_from_dendropy(self, tree: dendropy.Tree, method: str = 'pd') -> 'PDM':
        # Validate inputs.
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

        self._data = create_mat(len(self._indices), self._arr_default)

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

        return self

    @staticmethod
    def get_from_newick_file(path_tree: str, method: str = 'pd') -> 'PDM':
        tree = dendropy.Tree.get_from_path(path_tree,
                                           schema='newick',
                                           rooting='force-unrooted',
                                           preserve_underscores=True)
        return PDM.get_from_dendropy(tree, method)

    def save_to_path(self, path: str):
        with h5py.File(path, 'w') as f:
            f.create_dataset('indices',
                             data=[t.encode('ascii') for t in self._indices.get_keys()],
                             dtype=h5py.string_dtype(encoding='ascii'))
            f.create_dataset('data', data=self._data, chunks=True, dtype=self._d_type)
            f.create_dataset('tree_length', data=self._tree_length, dtype=self._d_type)
            f.create_dataset('arr_default', data=self._arr_default, dtype=self._d_type)

    @staticmethod
    def get_from_path(path: str) -> 'PDM':
        return PDM()._get_from_path(path)

    def _get_from_path(self, path: str) -> 'PDM':
        # TODO: STore the method and compare it
        self._indices = Indices()
        with h5py.File(path, 'r') as hf:
            self._arr_default = hf['arr_default'][()]
            self._data = hf['data'][()]
            self._tree_length = hf['tree_length'][()]
            for idx, key in enumerate(hf['indices'][()]):
                self._indices.add_key(key.decode('utf-8'))
        self._d_type = self._data.dtype
        return self

    def get_value(self, key_i: str, key_j: str, normalised: bool = False):
        value = self.get_value(key_i, key_j)
        if normalised:
            value *= 1 / self._tree_length
        return value
