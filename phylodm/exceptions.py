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


class PhyloDMException(Exception):
    """ Base exception for all PhyloDM exceptions thrown."""

    def __init__(self, message=''):
        Exception.__init__(self, message)


class DuplicateIndex(PhyloDMException):
    """Raised when a duplicate index is added to the Indices class."""

    def __init__(self, message=''):
        PhyloDMException.__init__(self, message)
