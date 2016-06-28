# gene_info

### Features
- Looks up canonical transcript properties for gene, using Ensembl API.
    * includes nucleotide and AA sequences, and CDS coordinate intervals.
- Coordinate conversion from hg19 position to CDS indices, and vice versa.
- Get nucleotide sequences for arbitrary hg19 intervals.
- Look up trinucleotide context (e.g. for Lego plot) and identify CpG sites
for mutation categorization and mutational signature quantification.
- Build BED files for silent/nonsilent transition/transversion sites.
    * These can be useful for generating custom coverage files for MutsigCV.

Compatible with Python 2 and 3.

## Usage

Look up gene using CanonicalInfo class.
```python
gene = gene_info.CanonicalInfo('TP53')
```

This will create an object with the following attributes:
```
hugo_symbol       'TP53'
chrom             '17'
gene_id           'ENSG00000141510'
transcript_id     'ENST00000269305'
strand            -1
n_exons           11
n_codons          394
cdna_len          2579
cds_len           1182
n_cds_intervals   10
cds_intervals     [(7579839, 7579912), (7579700, 7579721), ..., (7572927, 7573008)]
aa_seq            Seq('MEEPQSDPSVEPPLS...SD*', HasStopCodon(IUPACProtein(), '*'))
cds_seq           Seq('ATGGAGGAGCCGCAGTCAGATCCTAGCGT...TGA', IUPACUnambiguousDNA())
```

The initializer attempts to extract Ensembl transcript IDs and gene IDs from
 a bundled data file, with mappings from BioMart. Failed lookups (which may result
 from unofficial or out-of-date gene symbols) will fall back on the API.
  Alternatively, you can specify the ensembl transcript ID and/or gene ID if
  you have them.

### Coordinate mapping

Map from CDS index to genome coordinates:
```python
gene.get_hg_coord(0)
```
`>>> 7579912`

...and vice versa:
```python
gene.get_cds_index(7579912)
```
`>>> 0`

### Nucleotide sequence lookup

Provide a chromosome, start position and end position (both 1-based) for
nucleotide sequence lookup:
```python
gene_info.lookup_hg19('12', 25398285, 25398290)
```
```
>>> >12:25398285-25398290
    CAGCTC
```

### Mutation context

Get a basic categorization for a specified mutation as
{CpG, C:T, A:T} and {transition, transversion}, for use in MutSigCV, for example.
```python
gene_info.get_mutation_category('12', 25398285, 'C', 'A')
```
`>>> 'C:G_transversion'`

Alternatively, get the **trinucleotide context**, relative to the corresponding
pyrimidine mutation, used for [Lego plots](http://gdac.broadinstitute.org/runs/analyses__2016_01_28/reports/cancer/SKCM-TM/MutSigNozzleReportCV/SKCM-TM.final_analysis_set.maf.mutation_profile.png):
```python
gene_info.get_mutation_category_lego('12', 25398285, 'C', 'A')
```
`>>> ('C', 'C>A', 'A')`

### BED files

Obtain BED files for silent/nonsilent transition/transversion sites by providing an iterable of hugo symbols:
```python
gene_info.get_beds_for_hugo_list(['KRAS', 'TP53'])
```

This creates four BED files in the current directory, each containing
all corresponding sites for the specified genes:

- roi_ns_tv_final.bed
- roi_ns_ts_final.bed
- roi_s_ts_final.bed
- roi_s_tv_final.bed

You can specify paths for each file using method parameters.


## Installation

This package requires [bedtools](https://bedtools.readthedocs.io/en/latest/)
and a fasta file for hg19.

### Reference fasta

An hg19 reference fasta is available from the 1000 genomes project at:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/
- Download **human_g1k_v37.fasta{.gz, .fai}** to a directory of your choice
and unzip the .gz file.

### Configuring and installing

With bedtools installed and a fasta file downloaded, modify the `config.ini` file in
the directory *gene_info/config* to specify their paths.

Finally, from the root directory of the package, run `python setup.py install` to
install into your current Python environment.


## License

    Copyright (C) 2016 Stephen Gaffney

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
