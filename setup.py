from setuptools import setup

__author__ = 'Stephen G. Gaffney'

setup(name='gene_info',
      version='0.1',
      description='Fetch canonical transcript info and generate S/NS bed '
                  'files.',
      url='http://github.com/sggaffney',
      author='Stephen G. Gaffney',
      author_email='stephen.gaffney@yale.edu',
      license='GPLv3',
      packages=['gene_info'],
      package_data={'gene_info': ['data/*', 'config/*.ini']},
      install_requires=[
          'numpy', 'pyfaidx', 'requests', 'pandas', 'configparser', 'biopython'
      ],
      scripts=['bin/build_gene_beds', 'bin/lookup_hg19'],
      zip_safe=False,
      )
