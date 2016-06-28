#!/usr/bin/env python
import os
import subprocess

from gene_info import CanonicalInfo
from . import LookupFailedException, bedtools_dir


class TranscriptSiteInfo:
    """Populated by append_category_beds."""

    def __init__(self, hugo=None, chrom=None,
                 silent_ts=None, nonsilent_ts=None,
                 silent_tv=None, nonsilent_tv=None):
        self.hugo = hugo
        # self.refseq = refseq
        self.chrom = chrom
        self.site_list = list()
        self.path_dict = {
            ('transition', True): silent_ts,
            ('transition', False): nonsilent_ts,
            ('transversion', True): silent_tv,
            ('transversion', False): nonsilent_tv
            }

    def add_site(self, pos, kind=None, is_silent=True):
        """Add site category info to collection that will include whole
        transcript.
        :rtype : None
        """
        self.site_list.append((pos, kind, is_silent))

    def write_category_beds(self):
        """Write sites to file.
        file_dict = {
        ('transition', True): 'path_a',
        ('transition', False): 'path_b',
        ('transversion', True): 'path_c',
        ('transversion', False): 'path_d',
        }
        :rtype : None
        """
        streams = dict()  # 4 output streams
        try:
            for pair in self.path_dict:
                streams[pair] = open(self.path_dict[pair], 'a')

            for site in self.site_list:
                file_obj = streams[site[1:]]
                pos_str = "{}\t{}\t{}\t{}\n".\
                    format(self.chrom, site[0] - 1, site[0], self.hugo)
                file_obj.write(pos_str)
        finally:
            for stream in streams.values():
                stream.close()


def append_category_beds(hugo_symbol, silent_ts=None, silent_tv=None,
                         nonsilent_ts=None, nonsilent_tv=None):

    try:
        t = CanonicalInfo(hugo_symbol, lookup_seq=True)
        cds_seq = t.cds_seq
        aa_seq = t.aa_seq
        # g = GeneSeq(refseq_NM)
    except LookupFailedException as e:
        print(e.message)
        return
    # print(g.cds_seq)
    alt_dict = dict(A='CGT', C='AGT', G='ACT', T='ACG')
    transition_dict = {'A': 'G', 'C': 'T', 'G': 'A', 'T': 'C'}
    site_info = TranscriptSiteInfo(hugo=hugo_symbol,
                                   chrom=t.chrom,
                                   silent_ts=silent_ts,
                                   nonsilent_ts=nonsilent_ts,
                                   silent_tv=silent_tv,
                                   nonsilent_tv=nonsilent_tv)

    for codon_ind in range(len(aa_seq)):
        # identify previous and next base (prv_base, nxt_base)
        bases3 = cds_seq[codon_ind * 3:codon_ind * 3 + 3]
        aa_orig = aa_seq[codon_ind]
        # for each of 3 bases
        for subind in range(3):
            # has_silent = False
            # has_nonsilent = False
            cds_ind = codon_ind * 3 + subind
            hg_pos = t.get_hg_coord(cds_ind)
            # if bases3[subind].upper() in ['A', 'T']:
            # is_at = True
            #     is_cg = False
            #     is_cpg = False
            # else:
            #     is_at = False
            #     is_cpg = test_cpg(g.exonSet.chrom.replace('chr',''),hg_pos)
            #     assert isinstance(is_cpg, bool)
            #     is_cg = not is_cpg
            # for each alternative base
            this_base = bases3[subind]
            for alt in alt_dict[this_base]:  # 3 alternative bases
                if transition_dict[this_base].upper() == alt:
                    kind = 'transition'
                else:
                    kind = 'transversion'
                newbases = bases3.tomutable()
                newbases[subind] = alt
                new_aa = newbases.toseq().translate()
                if new_aa[0] == aa_orig:
                    site_info.add_site(hg_pos, kind=kind, is_silent=True)
                else:
                    site_info.add_site(hg_pos, kind=kind, is_silent=False)
    site_info.write_category_beds()
    # condense_beds_2cat(g.hugo)


def condense_beds_2cat(hugo):
    """Sort and merge both silent and nonsilent raw bed files for gene."""
    for categ in ['silent', 'nonsilent']:
        raw_path = "{}_{}.bed".format(hugo, categ)
        sorted_path = "{}_{}_sorted.bed".format(hugo, categ)
        condense_path = "{}_{}.final.bed".format(hugo, categ)

        cmd1 = ["sort", "-k1,1", "-k2,2n", raw_path, "-o", sorted_path]
        cmd2 = ["mergeBed", "-i", sorted_path, "-c", "4", "-o", "distinct"]
        subprocess.check_call(cmd1)
        with open(condense_path, 'w') as out:
            subprocess.check_call(cmd2, stdout=out)
        os.remove(raw_path)
        os.remove(sorted_path)


def condense_bed(raw_path):
    """Sort and merge both silent and nonsilent raw bed files for gene."""
    sorted_path = raw_path + "_sorted.bed"
    condense_path = raw_path + "_final.bed"

    cmd1 = ["sort", "-k1,1", "-k2,2n", raw_path, "-o", sorted_path]
    cmd2 = [os.path.join(bedtools_dir, "mergeBed"), "-i", sorted_path,
            "-c", "4", "-o", "distinct"]
    subprocess.check_call(cmd1)
    # with open(condense_path, 'w') as file:
    #     subprocess.check_call(cmd2, stdout=file)
    with open(os.devnull, "r") as fnullin:
        with open(condense_path, "w") as outfile:
            p = subprocess.Popen(cmd2, stdin=fnullin, stdout=outfile)
            p.wait()
    os.remove(raw_path)
    os.remove(sorted_path)


def get_beds_for_hugo_list(hugo_iterable,
                           silent_ts="roi_s_ts",
                           silent_tv="roi_s_tv",
                           nonsilent_ts="roi_ns_ts",
                           nonsilent_tv="roi_ns_tv"):
    """Takes list of refseq_NM strings and creates silent and nonsilent bed
    files."""
    path_list = [silent_ts, nonsilent_ts, silent_tv, nonsilent_tv]
    # for path in path_list:
    #     if os.path.exists(path):
    #         os.remove(path)

    for hugo in hugo_iterable:
        append_category_beds(hugo,
                             silent_ts=silent_ts,
                             nonsilent_ts=nonsilent_ts,
                             silent_tv=silent_tv,
                             nonsilent_tv=nonsilent_tv)
    for path in path_list:
        condense_bed(path)
