__author__ = 'sgg'

import unittest
from gene_info import CanonicalInfo

class CoordLookupTestCase(unittest.TestCase):
    """Tests for `gene_info_sgg.py`."""

    def test_start_cds(self):
        """Are coords correct for BRAF?"""
        t = CanonicalInfo('BRAF')
        self.assertTrue(t.get_cds_index(140624503) == 0)
        self.assertTrue(t.get_cds_index(140434397) == 2300)
        self.assertTrue(t.get_hg_coord(0) == 140624503)
        self.assertTrue(t.get_hg_coord(2300) == 140434397)

    def test_inconsistent_id(self):
        t = CanonicalInfo('MUC3A')
        self.assertTrue(t.gene_id == 'ENSG00000169894')
        self.assertTrue(t.transcript_id == 'ENST00000319509')
        self.assertTrue(t.cdna_len == 3338)

if __name__ == '__main__':
    unittest.main()
