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

import unittest

from phylodm.common import *

N_TESTS = 500


class TestCommon(unittest.TestCase):

    def test_create_mat_vector(self):
        default = 0
        for n in range(N_TESTS):
            vec = create_mat_vector(n, default)
            exp = len(np.triu_indices(n)[0])
            self.assertEqual(len(vec), exp)
            self.assertTrue(np.all(vec == default))

    def test_row_idx_from_mat_coords(self):
        expected = np.triu_indices(N_TESTS)
        for true_idx, (i, j) in enumerate(zip(expected[0], expected[1])):
            test_idx_a = row_idx_from_mat_coords(N_TESTS, i, j)
            test_idx_b = row_idx_from_mat_coords(N_TESTS, j, i)
            self.assertTrue(true_idx == test_idx_a == test_idx_b)

    def test_mat_shape_from_row_shape(self):
        for n in range(N_TESTS):
            row_len = len(np.triu_indices(n)[0])
            test = mat_shape_from_row_shape(row_len)
            self.assertTrue(n == test)

    def test_compact_int_mat(self):
        cases = (('int8', -128),
                 ('int16', -32768),
                 ('int32', -2147483648),
                 ('int64', -922337203685477),
                 ('uint8', 255),
                 ('uint16', 65535),
                 ('uint32', 4294967295),
                 ('uint64', 1844674407370955))
        for case, num in cases:
            true_arr = np.array([[1, 2], [3, num]], dtype=np.dtype('int64'))
            comp_arr = compact_int_mat(true_arr)
            self.assertEqual(case, comp_arr.dtype.name)
