from collections import deque

import dendropy
import h5py
from tqdm import tqdm

from phylodm.indices import Indices
from phylodm.symmat import SymMat


class PDM(SymMat):

    def __init__(self, n_indices: int, d_type):
        super().__init__(n_indices, d_type)

        self.tree_length = None

    def as_matrix(self, normalised=False):
        mat = super().as_matrix()
        if normalised:
            mat = mat * (1.0 / self.tree_length)
        return mat

    def process_dendropy_tree(self, tree: dendropy.Tree, method: str = 'pd'):
        # Validate inputs.
        if method == 'pd':
            use_pd = True
        elif method == 'node':
            use_pd = False
        else:
            raise NotImplemented(f'Unknown method: {method}')

        # Check that the end condition can be satisfied.
        if tree.seed_node.parent_node is not None:
            raise Exception('Dendropy did not seed the tree correctly.')

        # Store the tree length for normalising branch lengths.
        self.tree_length = tree.length()

        # Queue each of the leaf nodes for processing.
        queue = deque()
        for leaf_node in tree.leaf_node_iter():
            queue.append(leaf_node)

        # Create leaf indices
        for leaf in sorted(queue, key=lambda x: x.taxon.label):
            super()._indices.add_key(leaf.taxon.label)

        # Store lower triangle results in a vector.
        # n_leaf = len(self.leaf_idx)
        # d_type = np.float64 if use_pd else np.uint16
        # # self.results = np.zeros(int(n_leaf / 2 * (n_leaf - 1)), dtype=d_type)
        # self.results = np.zeros([n_leaf, n_leaf], dtype=d_type)

        # Track distances until the seed node, from each leaf node.
        parent_node = None
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
                    super().set_value(leaf_node.taxon.label, sister_node.taxon.label, sister_dist + tot_dist)
                    sisters_processed.add(sister_node)

                parent_node.child_dist[leaf_node] = tot_dist
                cur_node = parent_node

    def process_tree(self, path_tree, method: str = 'pd'):
        tree = dendropy.Tree.get_from_path(path_tree,
                                           schema='newick',
                                           rooting='force-unrooted',
                                           preserve_underscores=True)
        self.process_dendropy_tree(tree, method)

    def write(self, path: str):
        with h5py.File(path, 'w') as f:
            f.create_dataset('indices',
                             data=[t.encode('ascii') for t in self._indices.get_keys()],
                             dtype=h5py.string_dtype(encoding='ascii'))
            f.create_dataset('data', data=self._data, chunks=True)
            f.create_dataset('tree_length', data=self.tree_length)

    def read(self, path: str):
        self._indices = Indices()
        with h5py.File(path, 'r') as hf:
            self._data = hf['data'][()]
            self.tree_length = hf['tree_length'][()]
            for idx, key in enumerate(hf['indices'][()]):
                self._indices.add_key(key.decode('utf-8'))
