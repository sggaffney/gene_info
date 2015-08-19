#!/usr/bin/env python
import MySQLdb as mdb
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import subprocess
import os
import numpy as np
from ensembl_client import EnsemblRestClient

dbvars = {'host': 'localhost', 'db': 'CancerDB',
          'read_default_file': "~/.my.cnf"}

bed_bin_dir = os.path.join(os.path.expanduser('~'), 'bin', 'bed')  # bedtools


class NoIntervalsException(Exception):
    pass


class LookupFailedException(Exception):
    pass


class LengthMismatchException(Exception):
    pass


class CanonicalInfo:
    """Builds on ensembl api transcript info object."""
    def __init__(self, hugo_symbol):
        self.hugo_symbol = hugo_symbol
        transcript_gene = self.fetch_transcript_for_hugo(hugo_symbol)
        self.transcript_id, self.gene_id = transcript_gene
        self.cds_seq = self.get_cds_seq(self.transcript_id)
        self.aa_seq = self.cds_seq.translate()
        info = self.get_exons(self.transcript_id)
        self.chrom = info['seq_region_name']
        self.strand = info['strand']
        self.n_exons = len(info['Exon'])
        exon_starts = [i['start'] for i in info['Exon']]
        exon_ends = [i['end'] for i in info['Exon']]
        self.cdna_len = np.sum(np.array(exon_ends) - np.array(exon_starts) + 1)

        u5 = [i for i in info['UTR'] if i['object_type'] == 'five_prime_UTR']
        u3 = [i for i in info['UTR'] if i['object_type'] == 'three_prime_UTR']
        u5_starts = [i['start'] for i in u5]
        u5_ends = [i['end'] for i in u5]
        u3_starts = [i['start'] for i in u3]
        u3_ends = [i['end'] for i in u3]

        if self.strand == -1:
            u3_cutoff = max(u3_starts[:] + u3_ends[:])
            u5_cutoff = min(u5_starts[:] + u5_ends[:])
        else:
            u3_cutoff = min(u3_starts[:] + u3_ends[:])
            u5_cutoff = max(u5_starts[:] + u5_ends[:])

        cds_intervals = self._get_cds_intervals(self.strand, exon_starts,
                                                exon_ends, u5_cutoff, u3_cutoff)
        self.cds_intervals = cds_intervals
        self.n_cds_intervals = len(cds_intervals)
        a, b = zip(*cds_intervals)
        self.cds_len = sum(np.array(b) - np.array(a) + 1)
        if self.cds_len != len(self.cds_seq):
            raise LengthMismatchException("CDS intervals don't match seq"
                                          " length for {}".format(hugo_symbol))
        self.n_codons = self.cds_len/3

    def __repr__(self):
        return "<{classname} for {hugo}. {n_aa}aa, {cdslen}bp>".format(
            classname=self.__class__.__name__, hugo=self.hugo_symbol,
            n_aa=self.n_codons - 1, cdslen=self.cds_len)

    @staticmethod
    def get_exons(transcript_id, server='http://grch37.rest.ensembl.org'):
        """Lookup transcript, exon and UTR info from Ensembl Rest API."""
        client = EnsemblRestClient(server)
        transcript_info = client.perform_rest_action(
            '/lookup/id/{0}'.format(transcript_id),
            params={'object_type': 'transcript',
                    'expand': 1, 'utr': 1})
        return transcript_info

    @staticmethod
    def get_cds_seq(transcript_id, server='http://grch37.rest.ensembl.org'):
        """Lookup transcript, exon and UTR info from Ensembl Rest API."""
        client = EnsemblRestClient(server)
        seq_info = client.perform_rest_action(
            '/sequence/id/{0}'.format(transcript_id),
            params={'type': 'cds'})
        seq = Seq(seq_info['seq'], IUPAC.unambiguous_dna)
        return seq

    @staticmethod
    def fetch_transcript_for_hugo(hugo_symbol):
        """Get canonical transcript id and gene_id from ensembl."""
        con, rows, transcript_id, gene_id = None, None, None, None
        cmd1 = "SELECT `canonical_transcript`, `ensembl_gene_id` FROM " \
               "refs.`ensembl_canonical` WHERE hugo_symbol = {!r};".\
            format(hugo_symbol)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd1)
            if cur.rowcount != 1:
                raise LookupFailedException("Non-single row for transcript {}"
                                            .format(hugo_symbol))
            rows = cur.fetchall()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        # rows is [[transcript,gene],[transcript,gene],...]
        for row in rows:
            transcript_id = str(row[0])
            gene_id = str(row[1])
        return transcript_id, gene_id

    @staticmethod
    def _get_cds_intervals(strand, exon_starts, exon_ends,
                           u5_cutoff, u3_cutoff):
        """
        Get exon intervals after stripping out UTRs.
        for either strand, exons always listed in 5-3 order.
        start position always less than end position.
        """
        exon_intervals = zip(exon_starts, exon_ends)
        # force increasing coordinate exon order
        if strand == -1:
            exon_intervals.reverse()
            low_cutoff = u3_cutoff
            high_cutoff = u5_cutoff
        else:
            low_cutoff = u5_cutoff
            high_cutoff = u3_cutoff
        # remove low_cutoff section
        keep_intervals1 = []
        for interval in exon_intervals:
            if interval[0] > low_cutoff:
                keep_intervals1.append(interval)
                continue
            if interval[1] < low_cutoff:
                continue
            # now start <= cutoff <= end
            if low_cutoff+1 <= interval[1]:
                keep_intervals1.append((low_cutoff+1, interval[1]))
        keep_intervals2 = []
        # remove high_cutoff section
        for interval in keep_intervals1:
            if interval[1] < high_cutoff:
                keep_intervals2.append(interval)
                continue
            if interval[0] > high_cutoff:
                continue
            # now start <= cutoff <= end
            if high_cutoff-1 >= interval[0]:
                keep_intervals2.append((interval[0], high_cutoff-1))
        if strand == -1:
            keep_intervals2.reverse()
        return keep_intervals2

    def get_cds_index(self, hg_position):
        """Get coordinate, indexing from 0, in CDS, positive strand."""
        position = int(hg_position)  # ensure integer
        # get interval index containing position
        use_interval, interval_pos = None, None
        found_interval = False
        for ind, interval in enumerate(self.cds_intervals):
            if interval[0] <= position <= interval[1]:
                found_interval = True
                use_interval = ind
                if self.strand == 1:
                    interval_pos = position - interval[0]
                else:
                    interval_pos = interval[1] - position
                break
        if not found_interval:
            raise NoIntervalsException("No interval found for pos {} in {}".
                                       format(position, self.hugo_symbol))
        prev_length = sum([i[1] - i[0] + 1 for i in
                           self.cds_intervals[:use_interval]])
        # e.g. prev_length of 5 and interval_pos of 0 should give coord=5
        cds_index = interval_pos + prev_length
        return cds_index

    def get_hg_coord(self, cds_ind):
        """Get hg coordinate corresponding to cds_coordinate (0-indexed)."""
        cds_ind = int(cds_ind)  # ensure integer
        if cds_ind < 0:
            raise NoIntervalsException("No interval found for ind {} in {}".
                                       format(cds_ind, self.hugo_symbol))
        cds_lens = [i[1] - i[0] + 1 for i in self.cds_intervals]
        cds_lens.insert(0, 0)
        prv_cumul = np.cumsum(np.array(cds_lens))
        # position is in first interval where cumsum of previous intervals
        # is GREATER than cds_index
        try:
            interval_ind = np.nonzero(prv_cumul > cds_ind)[0][0] - 1
        except IndexError:
            raise NoIntervalsException("No interval found for ind {} in {}".
                                       format(cds_ind, self.hugo_symbol))
        subinterval_ind = cds_ind - prv_cumul[interval_ind]
        if self.strand == 1:
            hg_pos = self.cds_intervals[interval_ind][0] + subinterval_ind
        else:
            hg_pos = self.cds_intervals[interval_ind][1] - subinterval_ind
        return hg_pos


# import requests
# from cStringIO import StringIO
# from Bio import SeqIO
# def get_ensembl_seq(transcript_id):
#     """Return Biopython SeqIO sequence."""
#     params = {
#         'db': 'core',
#         'flank3_display': '0',
#         'flank5_display': '0',
#         'output': 'fasta',
#         'strand': 'feature',
#         't': transcript_id,
#         'param': 'cdna',
#         'genomic': 'off',
#         '_format': 'Text'
#         }
#     url = 'http://www.ensembl.org/Homo_sapiens/Export/Output/Transcript'
#     r = requests.get(url, params=params)
#     stream = StringIO(r.text)
#     return SeqIO.read(stream, "fasta")

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
        transcript."""
        self.site_list.append((pos, kind, is_silent))

    def write_category_beds(self):
        """Write sites to file.
        file_dict = {
        ('transition', True): 'path_a',
        ('transition', False): 'path_b',
        ('transversion', True): 'path_c',
        ('transversion', False): 'path_d',
        }
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
        t = CanonicalInfo(hugo_symbol)
        cds_seq = t.cds_seq
        aa_seq = t.aa_seq
        # g = GeneSeq(refseq_NM)
    except LookupFailedException as e:
        print e.message
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

    for codon_ind in xrange(len(aa_seq)):
        # identify previous and next base (prv_base, nxt_base)
        bases3 = cds_seq[codon_ind * 3:codon_ind * 3 + 3]
        aa_orig = aa_seq[codon_ind]
        # for each of 3 bases
        for subind in xrange(3):
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
        with open(condense_path, 'w') as file:
            subprocess.check_call(cmd2, stdout=file)
        os.remove(raw_path)
        os.remove(sorted_path)


def condense_bed(raw_path):
    """Sort and merge both silent and nonsilent raw bed files for gene."""
    sorted_path = raw_path + "_sorted.bed"
    condense_path = raw_path + "_final.bed"

    cmd1 = ["sort", "-k1,1", "-k2,2n", raw_path, "-o", sorted_path]
    cmd2 = [os.path.join(bed_bin_dir, "mergeBed"), "-i", sorted_path,
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


def get_beds_for_hugo_list(hugo_list,
                           silent_ts="/Users/Stephen/Downloads/roi_s_ts",
                           silent_tv="/Users/Stephen/Downloads/roi_s_tv",
                           nonsilent_ts="/Users/Stephen/Downloads/roi_ns_ts",
                           nonsilent_tv="/Users/Stephen/Downloads/roi_ns_tv"):
    """Takes list of refseq_NM strings and creates silent and nonsilent bed
    files."""
    path_list = [silent_ts, nonsilent_ts, silent_tv, nonsilent_tv]
    for path in path_list:
        if os.path.exists(path):
            os.remove(path)

    for hugo in hugo_list:
        append_category_beds(hugo,
                             silent_ts=silent_ts,
                             nonsilent_ts=nonsilent_ts,
                             silent_tv=silent_tv,
                             nonsilent_tv=nonsilent_tv)
    for path in path_list:
        condense_bed(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument("refseq", help="Refseq_NM.")
    parser.add_argument("hugo_file", help="File with list of hugo symbols.",
                        type=file)
    args = parser.parse_args()
    # refseq = args.refseq
    # append_category_beds(refseq)
    hugo_file = args.refseq_file
    hugos = list()
    try:
        for line in hugo_file:
            hugos.append(line.strip().split('\t')[1])
    finally:
        hugo_file.close()

    get_beds_for_hugo_list(hugos)

    # g = GeneSeq(refseq)
    # print g












