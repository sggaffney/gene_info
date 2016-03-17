#!/usr/bin/env python
"""
Requires /scratch/fasta/human_g1k_v37.fasta(.fai). pyfaidx for fast lookup.
"""

import os
import subprocess

from . import fasta_path
import pyfaidx

chroms = pyfaidx.Fasta(fasta_path, filt_function=lambda x: x[0] != 'G')


def lookup_hg19(chrom, start_pos, end_pos=None):
    """
    Returns bases for specified positions in hg19.
    If end_pos not provided, uses start_pos; i.e. returns one base.
    chrom should be of form 'X' not 'chrX'.
    """
    if not end_pos:
        end_pos = start_pos
    return str(chroms[chrom][start_pos-1:end_pos].seq)


def test_cpg(chrom, pos):
    """
    Return true if site is CpG. Tests adjacent sites looked up using hg19 fasta
    with samtools.
    """
    is_cpg = False  # set default is_cpg, will test and override.
    bases3 = str(chroms[chrom][pos-2:pos+1].seq)
    left = bases3[0]
    mid = bases3[1]
    right = bases3[2]
    if mid is 'C' and right is 'G':
        is_cpg = True
    elif mid is 'G' and left is 'C':
        is_cpg = True
    return is_cpg
