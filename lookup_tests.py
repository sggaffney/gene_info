__author__ = 'sgg'

import unittest
from gene_info_sgg import CanonicalInfo

class CoordLookupTestCase(unittest.TestCase):
    """Tests for `primes.py`."""

    def test_start_cds(self):
        """Is five successfully determined to be prime?"""
        t = CanonicalInfo('BRAF')
        self.assertTrue(t.get_cds_index(140624503) == 0)
        self.assertTrue(t.get_cds_index(140434397) == 2300)
        self.assertTrue(t.get_hg_coord(0) == 140624503)
        self.assertTrue(t.get_hg_coord(2300) == 140434397)

if __name__ == '__main__':
    unittest.main()
