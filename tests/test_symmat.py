import os
import shutil
import tempfile
import unittest

import numpy as np

from phylodm.symmat import SymMat


def generate_test_matrix(n):
    expected = np.full((n, n), -1)
    expected[np.triu_indices(n, 1)] = np.random.randint(0, 1e6, (n * (n - 1)) // 2)
    expected = expected + expected.T
    expected[np.diag_indices(n)] = np.random.randint(0, 1e6, n)
    assert (np.all(np.abs(expected - expected.T) < 1e-8))
    return expected


class TestSymMat(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='phylodm_test_')
        self.n = 500

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_sym_mat(self):
        mat = SymMat.get_from_shape(self.n, np.dtype('int32'), -1)
        expected = generate_test_matrix(self.n)

        for i in range(self.n):
            for j in range(i + 1):
                mat.set_value(str(i), str(j), expected[i][j])

        for i in range(self.n):
            for j in range(self.n):
                test_val = mat.get_value(str(i), str(j))
                true_val = expected[i][j]
                self.assertEqual(test_val, true_val)

    def test_as_matrix(self):
        sym_mat = SymMat.get_from_shape(self.n, np.dtype('int32'), -1)

        expected = generate_test_matrix(self.n)
        for i in range(self.n):
            for j in range(i + 1):
                sym_mat.set_value(str(i), str(j), expected[i][j])

        labels, mat = sym_mat.as_matrix()
        self.assertTrue(np.array_equal(expected, mat))

    def test_write(self):
        mat = SymMat.get_from_shape(self.n, np.dtype('int32'), -1)
        expected = generate_test_matrix(self.n)

        for i in range(self.n):
            for j in range(i + 1):
                mat.set_value(str(i), str(j), expected[i][j])

        # Write the matrix to disk.
        out_path = os.path.join(self.dir_tmp, 'output.h5py')
        mat.save_to_path(out_path)

        new_mat = SymMat.get_from_path(out_path)

        self.assertTrue(mat == new_mat)
