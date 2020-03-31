import os
import shutil
import tempfile
import unittest

import numpy as np

from phylodm.symmat import SymMat


class TestSymMat(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='phylodm_test_')
        self.path = os.path.join(self.dir_tmp, 'top_hit.tsv')

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_sym_mat(self):

        n = 1000
        mat = SymMat(n, np.int32, -1)

        expected = np.full((n, n), -1)
        expected[np.triu_indices(n, 1)] = np.random.randint(0, 1e6, (n * (n - 1)) // 2)
        expected = expected + expected.T
        expected[np.diag_indices(n)] = np.random.randint(0, 1e6, n)

        for i in range(n):
            for j in range(i + 1):
                mat.set_value(str(i), str(j), expected[i][j])

        for i in range(n):
            for j in range(n):
                test_val = mat.get_value(str(i), str(j))
                true_val = expected[i][j]
                self.assertEqual(test_val, true_val)

    def test_as_matrix(self):
        n = 1000
        sym_mat = SymMat(n, np.int32, -1)

        expected = np.full((n, n), -1)
        expected[np.triu_indices(n, 1)] = np.random.randint(0, 1e6, (n * (n - 1)) // 2)
        expected = expected + expected.T
        expected[np.diag_indices(n)] = np.random.randint(0, 1e6, n)

        for i in range(n):
            for j in range(i + 1):
                sym_mat.set_value(str(i), str(j), expected[i][j])

        mat = sym_mat.as_matrix()
        self.assertTrue(np.array_equal(expected, mat))

    def test_write(self):
        n = 1000
        mat = SymMat(n, np.int32, -1)

        expected = np.full((n, n), -1)
        expected[np.triu_indices(n, 1)] = np.random.randint(0, 1e6, (n * (n - 1)) // 2)
        expected = expected + expected.T
        expected[np.diag_indices(n)] = np.random.randint(0, 1e6, n)

        for i in range(n):
            for j in range(i + 1):
                mat.set_value(str(i), str(j), expected[i][j])

        # Write the matrix to disk.
        out = os.path.join(self.dir_tmp, 'output.h5py')
        mat.write(out)

        new_mat = SymMat(n, np.int32, -1)
        new_mat.read(out)

        self.assertTrue(mat == new_mat)
