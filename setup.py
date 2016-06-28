__author__ = 'Stephen G. Gaffney'

from setuptools import setup

setup(name='gene_info',
      version='0.1',
      description='Fetch canonical transcript info and generate S/NS bed files.',
      url='http://github.com/sggaffney',
      author='Stephen G. Gaffney',
      author_email='stephen.gaffney@yale.edu',
      license='GPLv3',
      packages=['gene_info'],
      install_requires=[
            'numpy', 'sqlalchemy', 'pyfaidx', 'requests', 'pandas',
            'configparser'
      ],
      scripts=['bin/build_gene_beds', 'bin/lookup_hg19'],
      zip_safe=True,
      # data_files=[('gene_matcher', ['gene_matcher/gene_lookup_refs.db'])]
      )
