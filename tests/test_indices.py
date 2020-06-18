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

from phylodm.exceptions import DuplicateIndex
from phylodm.indices import Indices


class TestIndices(unittest.TestCase):

    def setUp(self):
        self.indices = Indices()

    def tearDown(self):
        self.indices = None

    def test___init__(self):
        self.assertListEqual(self.indices._keys, [])
        self.assertDictEqual(self.indices._keys_idx, {})

    def test___eq__(self):
        # Basic conditions.
        a, b = Indices(), Indices()
        self.assertEqual(a, b)
        a.add_key('1')
        self.assertNotEqual(a, b)
        b.add_key('1')
        self.assertEqual(a, b)
        self.assertNotEqual(a, '1')

        # Check order equality.
        a.add_keys(['2', '3'])
        b.add_keys(['3', '2'])
        self.assertNotEqual(a, b)

    def test___len__(self):
        self.assertEqual(len(self.indices), 0)
        self.indices.add_key('1')
        self.assertEqual(len(self.indices), 1)
        self.indices.add_keys(['2', '3', '4', '5', '6'])
        self.assertEqual(len(self.indices), 6)

    def test___contains__(self):
        self.assertFalse('foo' in self.indices)
        self.indices.add_key('foo')
        self.assertTrue('foo' in self.indices)

    def test_get_keys(self):
        self.assertTupleEqual(tuple([]), self.indices.get_keys())
        self.indices._keys = ['foo', 'bar']
        self.indices._keys_idx = {'foo': 0, 'bar': 1}
        self.assertTupleEqual(tuple(['foo', 'bar']), self.indices.get_keys())

    def test_add_key(self):
        idx_foo = self.indices.add_key('foo')
        idx_bar = self.indices.add_key('bar')
        self.assertTupleEqual(tuple([idx_foo, idx_bar]), tuple([0, 1]))
        self.assertTupleEqual(tuple(['foo', 'bar']), self.indices.get_keys())
        self.assertDictEqual({'foo': 0, 'bar': 1}, self.indices.get_key_indices())
        self.assertRaises(DuplicateIndex, self.indices.add_key, 'foo')

    def test_add_keys(self):
        idx = self.indices.add_keys(['foo', 'bar'])
        self.assertTupleEqual(tuple(idx), tuple([0, 1]))
        self.assertTupleEqual(tuple(['foo', 'bar']), self.indices.get_keys())
        self.assertDictEqual({'foo': 0, 'bar': 1}, self.indices.get_key_indices())

    def test_get_key_indices(self):
        self.assertTupleEqual(tuple([]), self.indices.get_keys())
        self.indices._keys = ['foo', 'bar']
        self.indices._keys_idx = {'foo': 0, 'bar': 1}
        self.assertDictEqual({'foo': 0, 'bar': 1}, self.indices.get_key_indices())

    def test_get_key_idx(self):
        self.assertTupleEqual(tuple([]), self.indices.get_keys())
        self.indices._keys = ['foo', 'bar']
        self.indices._keys_idx = {'foo': 0, 'bar': 1}
        self.assertEqual(0, self.indices.get_key_idx('foo'))
        self.assertEqual(1, self.indices.get_key_idx('bar'))
