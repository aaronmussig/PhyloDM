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

import os
import shutil
import tempfile
import unittest

import numpy as np

from phylodm.symmat import SymMat

N_TESTS = 100


def generate_test_matrix(n):
    expected = np.full((n, n), -1, dtype=np.int32)
    expected[np.triu_indices(n, 1)] = np.random.randint(0, 1e6, (n * (n - 1)) // 2)
    expected = expected + expected.T
    expected[np.diag_indices(n)] = np.random.randint(0, 1e6, n)
    assert (np.all(np.abs(expected - expected.T) < 1e-8))
    return expected


class TestSymMat(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='phylodm_test_')
        self.n = N_TESTS

        # Create a test matrix and populate it with expected data.
        self.expected = generate_test_matrix(self.n)
        self.mat = SymMat.get_from_indices(list(map(str, range(self.n))), np.dtype('int32'), -1)
        for i in range(self.n):
            for j in range(i + 1):
                self.mat.set_value(str(i), str(j), self.expected[i][j])

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test___eq__(self):
        self.assertEqual(self.mat, self.mat)
        self.assertNotEqual(self.mat, 2)
        self.assertNotEqual(self.mat, SymMat.get_from_shape(self.n, np.dtype('int32'), -1))

    def test_sym_mat(self):
        for i in range(self.n):
            for j in range(self.n):
                test_val = self.mat.get_value(str(i), str(j))
                true_val = self.expected[i][j]
                self.assertEqual(test_val, true_val)

    def test_as_matrix(self):
        labels, mat = self.mat.as_matrix()
        self.assertTrue(np.array_equal(self.expected, mat))

    def test_write(self):
        out_path = os.path.join(self.dir_tmp, 'output.mat')
        self.mat.save_to_path(out_path)

        new_mat = SymMat.get_from_path(out_path)
        self.assertTrue(self.mat == new_mat)

    def test_remove_keys(self):

        # Generate a bunch of test conditions to simulate, with various keys missing.
        min_dim, max_dim = 3, 30
        for n in range(min_dim, max_dim):
            # Remove until there are only two dimensions left.
            for m in range(n - 2):

                # Determine what indices will be removed the test matrix.
                start_indices = set(range(n))
                remove_indices = {int(x) for x in np.random.permutation(list(start_indices))[0:m + 1]}
                idx_after_remove = start_indices - remove_indices
                true_vals = generate_test_matrix(n)

                # Create the full-sized test matrix.
                test_mat = SymMat.get_from_indices(sorted(list(map(str, start_indices))), np.dtype('int32'), -1)
                true_mat = SymMat.get_from_indices(sorted(list(map(str, idx_after_remove))), np.dtype('int32'), -1)
                for i in range(n):
                    for j in range(i + 1):
                        test_mat.set_value(str(i), str(j), true_vals[i][j])
                        if i not in remove_indices and j not in remove_indices:
                            true_mat.set_value(str(i), str(j), true_vals[i][j])

                # Remove the indices
                test_mat.remove_keys(sorted([str(x) for x in remove_indices]))
                self.assertEqual(test_mat, true_mat)
