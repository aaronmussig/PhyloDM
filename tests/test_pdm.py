import os
import shutil
import tempfile
import unittest

import numpy as np
from dendropy.simulate import treesim

from phylodm.pdm import PDM


def get_test_tree(n: int) -> dict:
    tree = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=n)
    pdm = tree.phylogenetic_distance_matrix()
    taxa = sorted(pdm.taxon_iter())
    pd_mat = np.zeros((n, n))
    pd_mat_norm = np.zeros((n, n))
    nd_mat = np.zeros((n, n))
    nd_mat_norm = np.zeros((n, n))
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
            'taxa': tuple([x.label for x in taxa]),
            'pd_mat': pd_mat,
            'pd_mat_norm': pd_mat_norm,
            'nd_mat': nd_mat,
            'nd_mat_norm': nd_mat_norm}


class TestPDM(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='phylodm_test_')
        self.true_data = get_test_tree(500)

        # Write the file for testing the newick get method.
        self.path_newick = os.path.join(self.dir_tmp, 'newick.tree')
        self.true_data['tree'].write(path=self.path_newick, schema='newick')

        # Save each of the methods to disk.
        self.path_pd = os.path.join(self.dir_tmp, 'pd.mat')
        self.path_nd = os.path.join(self.dir_tmp, 'nd.mat')
        PDM.get_from_dendropy(self.true_data['tree'], method='pd').save_to_path(self.path_pd)
        PDM.get_from_dendropy(self.true_data['tree'], method='node').save_to_path(self.path_nd)

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_pd_get_methods(self):
        pdm_methods = {'get_from_dendropy': PDM.get_from_dendropy(self.true_data['tree'], method='pd'),
                       'get_from_newick_file': PDM.get_from_newick_file(self.path_newick, method='pd'),
                       'get_from_path': PDM.get_from_path(self.path_pd)}

        for method, pdm in pdm_methods.items():
            for normalised in (False, True):

                labels, mat = pdm.as_matrix(normalised=normalised)
                self.assertAlmostEqual(self.true_data['tree'].length(), pdm._tree_length, 10)
                self.assertEqual(pdm._d_type, np.dtype('float64'))
                self.assertEqual(pdm._arr_default, 0.0)
                self.assertSetEqual(set(labels), set(self.true_data['taxa']))

                if normalised:
                    self.assertTrue(np.all(np.abs(mat - self.true_data['pd_mat_norm']) < 1e-8))
                else:
                    self.assertTrue(np.all(np.abs(mat - self.true_data['pd_mat']) < 1e-8))

    def test_node_get_methods(self):
        pdm_methods = {'get_from_dendropy': PDM.get_from_dendropy(self.true_data['tree'], method='node'),
                       'get_from_newick_file': PDM.get_from_newick_file(self.path_newick, method='node'),
                       'get_from_path': PDM.get_from_path(self.path_nd)}

        for method, pdm in pdm_methods.items():
            for normalised in (False, True):

                labels, mat = pdm.as_matrix(normalised=normalised)
                self.assertEqual(len(self.true_data['tree'].edges()), pdm._tree_length)
                self.assertEqual(pdm._d_type, np.dtype('uint32'))
                self.assertEqual(pdm._arr_default, 0)
                self.assertSetEqual(set(labels), set(self.true_data['taxa']))

                if normalised:
                    self.assertTrue(np.all(np.abs(mat - self.true_data['nd_mat_norm']) < 1e-8))
                else:
                    self.assertTrue(np.all(np.abs(mat - self.true_data['nd_mat']) < 1e-8))
