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

from typing import Dict, List


class Indices(object):

    def __init__(self):
        self._keys = list()
        self._keys_idx = dict()

    def __eq__(self, other):
        if isinstance(other, Indices):
            return self._keys == other._keys and \
                   self._keys_idx == other._keys_idx
        return False

    def get_keys(self) -> tuple:
        return tuple(self._keys)

    def add_key(self, key: str) -> int:
        if key in self._keys_idx:
            raise ValueError(f'Cannot add duplicate key: {key}')
        new_idx = len(self._keys)
        self._keys.append(key)
        self._keys_idx[key] = new_idx
        return new_idx

    def add_keys(self, keys: List[str]):
        for key in keys:
            self.add_key(key)

    def get_key_indices(self) -> Dict[str, int]:
        return self._keys_idx

    def get_key_idx(self, key: str) -> int:
        return self._keys_idx[key]

    def contains(self, key: str) -> bool:
        return key in self._keys_idx
