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

from typing import Dict, Collection, Tuple, List

from phylodm.exceptions import DuplicateIndex


class Indices(object):
    """Stores a unique list of strings, preserving their order and providing
    constant time lookup for the order in which they were added."""

    def __init__(self):
        """Instance variables are populated via the add key(s) method."""
        self._keys: List[str] = list()
        self._keys_idx: Dict[str, int] = dict()

    def __eq__(self, other) -> bool:
        """Equal if the indices and the positions are the same."""
        if isinstance(other, Indices):
            return self._keys == other._keys and \
                   self._keys_idx == other._keys_idx
        return False

    def __len__(self) -> int:
        """The number of indices stored."""
        return len(self._keys)

    def __contains__(self, key: str) -> bool:
        """Returns True if a key is in the list of indices."""
        return key in self._keys_idx

    def add_key(self, key: str):
        """Add a key, a DuplicateIndex is raised if the key is not unique."""
        if key in self._keys_idx:
            raise DuplicateIndex(f'Cannot add duplicate key: {key}')
        new_idx = len(self._keys)
        self._keys.append(key)
        self._keys_idx[key] = new_idx
        return new_idx

    def add_keys(self, keys: Collection[str]) -> List[int]:
        """Adds a collection of keys to the indices."""
        indexes = list()
        for key in keys:
            indexes.append(self.add_key(key))
        return indexes

    def get_keys(self) -> Tuple[str]:
        """Return a tuple of the indices stored."""
        return tuple(self._keys)

    def get_key_idx(self, key: str) -> int:
        """Returns the index of the specified key."""
        return self._keys_idx[key]

    def get_key_indices(self) -> Dict[str, int]:
        """Returns a dictionary of all keys and their indices."""
        return self._keys_idx
