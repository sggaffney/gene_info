__author__ = 'Stephen G. Gaffney'

"""
DOWNLOAD human_g1k_v37.fasta and human_g1k_v37.fasta.fai from:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/.
(SEE http://www.1000genomes.org/category/assembly for description.)

MOVE TO /scratch/fasta, or change path below.
"""

import os

import sqlalchemy


# CREATE SQLALCHEMY DB ENGINE
home_dir = os.path.expanduser('~')
cnf_path = os.path.join(home_dir, '.my.cnf')
db_url = sqlalchemy.engine.url.URL(drivername='mysql', host='localhost',
             database='refs',
             query={'read_default_file': cnf_path})
engine = sqlalchemy.create_engine(name_or_url=db_url)

fasta_path = '/scratch/fasta/human_g1k_v37.fasta'


class NoIntervalsException(Exception):
    pass


class LookupFailedException(Exception):
    pass

from gene_info import CanonicalInfo
from lookup_hg19 import lookup_hg19, test_cpg
from categorize import get_mutation_category
