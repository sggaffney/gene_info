import os
import sys
import pandas as pd

from configparser import ConfigParser

"""
DOWNLOAD human_g1k_v37.fasta and human_g1k_v37.fasta.fai from:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/.
(SEE http://www.1000genomes.org/category/assembly for description.)

SPECIFY path (hg19_fasta_path) in config.ini file.
"""

__author__ = 'Stephen G. Gaffney'

# get fasta path and bedtools dir from config.ini
parser = ConfigParser()
config_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'config', 'config.ini')
parser.read(config_path)

fasta_path = parser.get('default', 'hg19_fasta_path')
if not os.path.isfile(fasta_path):
    sys.exit("Please update hg19_fasta_path in {}".format(config_path))

bedtools_dir = parser.get('default', 'bedtools_dir')
if not os.path.isdir(bedtools_dir):
    sys.exit("Please update bedtools_dir in {}".format(config_path))

ensembl_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'data', 'ensembl_canonical.txt')
ensembl_df = pd.read_csv(ensembl_path, sep='\t', index_col=0)


class NoIntervalsException(Exception):
    pass


class LookupFailedException(Exception):
    pass


from gene_info import CanonicalInfo
from lookup_hg19 import lookup_hg19, test_cpg
from categorize import get_mutation_category, get_mutation_category_lego
from build_gene_beds import get_beds_for_hugo_list
