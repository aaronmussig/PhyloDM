from __future__ import annotations

from typing import Optional, List
from typing import TYPE_CHECKING

import dendropy

if TYPE_CHECKING:
    import numpy as np

from .pdm import PhyloDM as PDM


class PhyloDM:

    def __init__(self):
        """Initialize the PhyloDM object."""
        self._rs = PDM()

    @classmethod
    def load_from_newick_path(cls, path: str) -> 'PhyloDM':
        """Load a tree from a Newick file.

        Args:
            path: The path to the Newick file.
        """
        try:
            pdm = cls()
            pdm._rs.load_from_newick_path(path=path)
            return pdm
        except Exception as e:
            print(f'Unable to load newick tree using light_phylogeny (Rust). '
                  f'This is likely due to it not supporting the extended Newick format... '
                  f'falling back to DendroPy to load the tree.')
            tree = dendropy.Tree.get(path=path, schema='newick')
            return cls.load_from_dendropy(tree)

    @classmethod
    def load_from_dendropy(cls, tree: dendropy.Tree) -> 'PhyloDM':
        """Load a tree from a Dendropy tree object.

        Args:
            tree: The Dendropy tree object.
        """
        pdm = cls()

        node_to_id = dict()
        for node in tree.postorder_node_iter():
            if node.taxon and node.taxon.label:
                new_node_id = pdm.add_node(taxon=node.taxon.label)
            else:
                new_node_id = pdm.add_node()
            node_to_id[node] = new_node_id

        for node in tree.postorder_node_iter():
            if node.parent_node is not None:
                pdm.add_edge(parent_id=node_to_id[node.parent_node],
                             child_id=node_to_id[node],
                             length=node.edge_length)
        return pdm

    def add_node(self, taxon: Optional[str] = None) -> int:
        """Add a new node to the tree.

        Args:
            taxon: The taxon name if this is a leaf node.

        Returns:
            The index of the new node.
        """
        return self._rs.add_node(taxon=taxon)

    def add_edge(self, parent_id: int, child_id: int, length: float):
        """Add an edge between the two nodes

        Args:
            parent_id: The index of the parent node.
            child_id: The index of the child node.
            length: The length of the edge.
        """
        return self._rs.add_edge(parent_id=parent_id, child_id=child_id, length=length)

    def get_nodes(self) -> List[int]:
        """Return all node indexes in the tree."""
        return self._rs.get_nodes()

    def dm(self, norm: Optional[bool] = False) -> np.ndarray:
        """Returns a symmetrical distance matrix.

        Args:
            norm: If True, the matrix is normalized by branch length.
        """
        return self._rs.dm(norm=norm)

    def taxa(self) -> List[str]:
        """Returns a list of all taxa within the tree."""
        return self._rs.taxa()

    def length(self) -> float:
        """Returns the total length of the tree (sum of branch lengths)."""
        return self._rs.length()

    def compute_row_vec(self):
        """Compute the row vector for the tree (required if not initialised from a Newick file)."""
        return self._rs.compute_row_vec()

    def distance(self, a: str, b: str, norm: Optional[bool] = False) -> float:
        """Compute the distance between two taxa.

        Args:
            a: The first taxon.
            b: The second taxon.
            norm: If the distance should be normalised by the sum of branch lengths.

        Returns:
            The distance between the two taxa.
        """
        return self._rs.distance(a=a, b=b, norm=norm)
