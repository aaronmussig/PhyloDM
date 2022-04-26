import os
import tempfile
import unittest

import dendropy
import numpy as np
from dendropy.simulate import treesim

from phylodm import PhyloDM


def add_trifurication(tree):
    parent_node = list(tree.leaf_node_iter())[0].parent_node

    t1 = dendropy.Taxon(f'X1')
    t2 = dendropy.Taxon(f'X2')
    t3 = dendropy.Taxon(f'X3')

    tree.taxon_namespace.add_taxon(t1)
    tree.taxon_namespace.add_taxon(t2)
    tree.taxon_namespace.add_taxon(t3)

    child_a = dendropy.Node(edge_length=1.234)
    child_b = dendropy.Node(edge_length=1.234)
    child_c = dendropy.Node(edge_length=4.123)

    child_a.taxon = t1
    child_b.taxon = t2
    child_c.taxon = t3

    parent_node.add_child(child_a)
    parent_node.add_child(child_b)
    parent_node.add_child(child_c)


def get_test_tree(n: int, trifurication=False) -> dict:
    tree = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, num_extant_tips=n)

    if trifurication:
        add_trifurication(tree)

    for i, edge in enumerate(tree.preorder_node_iter()):
        edge.edge_length = i

    n_taxa = len(tree.taxon_namespace)
    pdm = tree.phylogenetic_distance_matrix()
    taxa = sorted(pdm.taxon_iter())
    pd_mat = np.zeros((n_taxa, n_taxa))
    pd_mat_norm = np.zeros((n_taxa, n_taxa))
    nd_mat = np.zeros((n_taxa, n_taxa))
    nd_mat_norm = np.zeros((n_taxa, n_taxa))
    for i, t1 in enumerate(taxa):
        for j, t2 in enumerate(taxa):
            pd = pdm.patristic_distance(t1, t2)
            pd_norm = pdm.patristic_distance(t1, t2, is_normalize_by_tree_size=True)
            nd = pdm.path_edge_count(t1, t2)
            nd_norm = pdm.path_edge_count(t1, t2, is_normalize_by_tree_size=True)
            pd_mat[i, j] = pd
            pd_mat_norm[i, j] = pd_norm
            nd_mat[i, j] = nd
            nd_mat_norm[i, j] = nd_norm

    return {'tree': tree,
            'length': tree.length(),
            'taxa': tuple([x.label for x in taxa]),
            'pd_mat': pd_mat,
            'pd_mat_norm': pd_mat_norm,
            'nd_mat': nd_mat,
            'nd_mat_norm': nd_mat_norm}


class TestPhyloDM(unittest.TestCase):

    def test_add_node(self):
        pdm = PhyloDM()
        pdm.add_node(taxon=None)
        pdm.add_node(taxon='b')

        nodes = pdm.get_nodes()
        self.assertEqual(len(nodes), 2)

    def test_add_edge(self):
        pdm = PhyloDM()
        a = pdm.add_node(taxon=None)
        b = pdm.add_node(taxon='b')
        pdm.add_edge(a, b, length=0.2)

    def test_dm(self):
        test_tree = get_test_tree(50)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = os.path.join(tmpdir, 'test.tree')
            with open(tmp_path, 'w') as f:
                f.write(test_tree['tree'].as_string(schema='newick')[5:])

            pdm = PhyloDM.load_from_newick_path(tmp_path)

        dm = pdm.dm(norm=False)
        delta = test_tree['pd_mat'] - dm

        self.assertTrue(np.all(delta < 1e-6))
        self.assertAlmostEqual(pdm.length(), test_tree['length'], places=6)
        self.assertTrue(test_tree['taxa'] == tuple(pdm.taxa()))

    def test_dm_norm(self):
        test_tree = get_test_tree(50)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = os.path.join(tmpdir, 'test.tree')
            with open(tmp_path, 'w') as f:
                f.write(test_tree['tree'].as_string(schema='newick')[5:])

            pdm = PhyloDM.load_from_newick_path(tmp_path)

        dm_norm = pdm.dm(norm=True)
        delta_norm = test_tree['pd_mat'] / pdm.length() - dm_norm

        self.assertTrue(np.all(delta_norm < 1e-6))
        self.assertAlmostEqual(pdm.length(), test_tree['length'], places=6)
        self.assertTrue(test_tree['taxa'] == tuple(pdm.taxa()))

    def test_load_from_dendropy(self):
        test_tree = get_test_tree(50)

        pdm = PhyloDM.load_from_dendropy(test_tree['tree'])
        dm = pdm.dm(norm=False)
        delta = test_tree['pd_mat'] - dm

        self.assertTrue(np.all(delta < 1e-6))
        self.assertAlmostEqual(pdm.length(), test_tree['length'], places=6)
        return

    def test_load_from_dendropy_with_trifurication(self):
        test_tree = get_test_tree(50, trifurication=True)

        pdm = PhyloDM.load_from_dendropy(test_tree['tree'])
        dm = pdm.dm(norm=False)
        delta = test_tree['pd_mat'] - dm

        self.assertTrue(np.all(delta < 1e-6))
        self.assertAlmostEqual(pdm.length(), test_tree['length'], places=6)
        self.assertTrue(test_tree['taxa'] == tuple(pdm.taxa()))
        return
