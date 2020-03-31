import unittest

from phylodm.indices import Indices


class TestIndices(unittest.TestCase):

    def setUp(self):
        self.indices = Indices()

    def tearDown(self):
        self.indices = None

    def test___init__(self):
        self.assertListEqual(self.indices._keys, [])
        self.assertDictEqual(self.indices._keys_idx, {})

    def test_get_keys(self):
        self.assertTupleEqual(tuple([]), self.indices.get_keys())
        self.indices._keys = ['foo', 'bar']
        self.indices._keys_idx = {'foo': 0, 'bar': 1}
        self.assertTupleEqual(tuple(['foo', 'bar']), self.indices.get_keys())

    def test_add_key(self):
        self.indices.add_key('foo')
        self.indices.add_key('bar')
        self.assertTupleEqual(tuple(['foo', 'bar']), self.indices.get_keys())
        self.assertDictEqual({'foo': 0, 'bar': 1}, self.indices.get_key_indices())
        self.assertRaises(ValueError, self.indices.add_key, 'foo')

    def test_add_keys(self):
        self.indices.add_keys(['foo', 'bar'])
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

    def test_contains(self):
        self.assertFalse(self.indices.contains('foo'))
        self.indices.add_key('foo')
        self.assertTrue(self.indices.contains('foo'))
