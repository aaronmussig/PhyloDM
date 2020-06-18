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
    expected = np.full((n, n), -1)
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
