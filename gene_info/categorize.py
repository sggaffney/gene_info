from gene_info.lookup_hg19 import test_cpg, get_seq_triplet

__author__ = 'sgg'


class BadRefException(Exception):
    pass


transition_dict = {'A': 'G', 'C': 'T', 'G': 'A', 'T': 'C'}

change_map = {'G>A': 'C>T',
              'G>T': 'C>A',
              'G>C': 'C>G',
              'T>C': 'A>G',
              'T>G': 'A>C',
              'T>A': 'A>T'}

lego_order = ['TT', 'CT', 'AT', 'GT', 'TC', 'CC', 'AC', 'GC', 'TA', 'CA', 'AA',
              'GA', 'TG', 'CG', 'AG', 'GG']


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
    except KeyError as e:
        print(e.message)
        return 'unknown'

    base3 = get_seq_triplet(chrom, pos)
    true_ref = str(base3)[1]
    if ref != true_ref:
        print('chr{chrom}:{pos} ref is {actual} not {ref}'.
              format(chrom=chrom, pos=pos, actual=true_ref, ref=ref))

    if transition_dict[ref] == alt:
        change = 'transition'
    else:
        change = 'transversion'

    if ref in ['C', 'G']:
        site = 'CpG' if is_cpg else 'C:G'
    else:
        site = 'A:T'

    return '_'.join([site, change])


def get_mutation_category_lego(chrom, pos, ref, alt):
    """Get base change (6 categories), base_before and base_after.

    Args: chrom must be str: 1-22,X,Y,MT,G...

    Returns: (tuple) base_before, change_str, base_after.
        e.g. ('C', 'C>A', 'A') for KRAS G12C mutation
    """
    if len(ref) > 1 and len(ref) == len(alt):
        return None, 'multiple', None
    if len(ref) != 1 or len(alt) != 1 or '-' in [ref, alt]:
        return None, 'indel', None

    base3 = get_seq_triplet(chrom, pos)
    change_str = '>'.join([ref, alt])
    # get reverse complement if dealing with 6 'flipped' mutations
    true_ref = str(base3)[1]
    if ref != true_ref:
        print('chr{chrom}:{pos} ref is {actual} not {ref}'.
              format(chrom=chrom, pos=pos, actual=true_ref, ref=ref))
    if change_str in change_map:
        change_str = change_map[change_str]
        base3 = -base3
    base3 = str(base3)

    base_before = base3[0]
    base_after = base3[2]
    return base_before, change_str, base_after
