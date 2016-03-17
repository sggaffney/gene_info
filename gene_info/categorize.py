__author__ = 'sgg'

from lookup_hg19 import test_cpg


transition_dict = {'A': 'G', 'C': 'T', 'G': 'A', 'T': 'C'}


def get_mutation_category(chrom, pos, ref, alt):
    """Get mutation category: A:T, C:G, Cpg transition or transversion.

    Returns (str): one of 'A:T_transition', 'A:T_transversion', ...
    """
    chrom, ref, alt = [i.upper() for i in [chrom, ref, alt]]
    if len(ref) > 1 and len(ref) == len(alt):
        return 'multiple'
    if len(ref) != 1 or len(alt) != 1 or '-' in [ref, alt]:
        return 'indel'
    try:
        is_cpg = test_cpg(chrom, pos)
    except KeyError, e:
        print e.message
        return 'unknown'

    if transition_dict[ref] == alt:
        change = 'transition'
    else:
        change = 'transversion'

    if ref in ['C', 'G']:
        site = 'CpG' if is_cpg else 'C:G'
    else:
        site = 'A:T'

    return '_'.join([site, change])
