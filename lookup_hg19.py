#!/usr/bin/env python
""" Requires ~/ref/human_g1k_v37.fasta(.fai) and samtools for lookup using faidx
"""

import subprocess
import os
import sys

def lookup_hg19(chrom, start_pos, end_pos=None):
    """
    Returns bases for specified positions in hg19.
    If end_pos not provided, uses start_pos; i.e. returns one base.
    chrom should be of form 'X' not 'chrX'.
    """
    if not end_pos:
        end_pos = start_pos
    user_dir = os.path.expanduser('~')
    fasta_path = os.path.join(user_dir, 'ref', 'human_g1k_v37.fasta')
    samtools_path = os.path.join(user_dir,'bin','samtools_bin','samtools')
    cmd = "{samtools} faidx {fastapath} {chrom}:{start_pos}-{end_pos}".format(
            samtools=samtools_path, fastapath=fasta_path, chrom=chrom,
            start_pos=start_pos, end_pos=end_pos)
    str_out = subprocess.check_output(cmd.split(' '))
    out_lines = str_out.split('\n')
    bases = ''.join(out_lines[1:])
    return bases

def test_cpg(chrom, pos):
    """
    Return true if site is CpG. Tests adjacent sites looked up using hg19 fasta with samtools.
    """
    is_cpg = False  # set default is_cpg, will test and override.
    user_dir = os.path.expanduser('~')
    fasta_path = os.path.join(user_dir, 'ref', 'human_g1k_v37.fasta')
    samtools_path = os.path.join(user_dir,'bin','samtools_bin','samtools')
    cmd = "{samtools} faidx {fastapath} {chrom}:{posA}-{posB}".format(
            samtools=samtools_path, fastapath=fasta_path, chrom=chrom,
            posA=int(pos)-1, posB=int(pos)+1)
    str_out = subprocess.check_output(cmd.split(' '))
    bases3 = str_out.split('\n')[1]
    left = bases3[0]
    mid = bases3[1]
    right = bases3[2]
    if mid is 'C' and right is 'G':
        is_cpg = True
    elif mid is 'G' and left is 'C':
        is_cpg = True
    return is_cpg

if __name__ == '__main__':
    if len(sys.argv)>1:
        chrom = sys.argv[1]
        pos = sys.argv[2]
    else:
        chrom, pos = '12', 25380289
    bases_str = lookup_hg19(chrom, pos)
    print(bases_str)
    print('is_cpg = {}'.format(test_cpg(chrom,pos)))