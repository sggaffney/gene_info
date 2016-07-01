"""Uses Ensembl API to look up canonical transcript info for given gene."""

import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from gene_info.ensembl_client import EnsemblRestClient, LookupFailedException
from . import NoIntervalsException, ensembl_df


client = EnsemblRestClient(server='http://grch37.rest.ensembl.org')


class LengthMismatchException(Exception):
    pass


class CanonicalInfo(object):
    """Builds on ensembl api transcript info object.

    Provide hugo_symbol and (optionally) ensembl transcript and gene to cut down
    on lookups. Specify lookup_seq=True for cds and aa sequence lookups.
    """
    def __init__(self, hugo_symbol=None, transcript_id=None, ensembl_gene=None,
                 lookup_seq=False):
        """"""
        self.hugo_symbol = hugo_symbol
        self.transcript_id = transcript_id
        self.gene_id = ensembl_gene
        self._cds_seq = None
        self._aa_seq = None
        info = None

        try:
            # fetch_ids_for_hugo OR fetch_transcript could fail
            if transcript_id is None:
                transcript_gene = fetch_ids_for_hugo(hugo_symbol)
                self.transcript_id, self.gene_id = transcript_gene
            if self.transcript_id is not None:
                # REST LOOKUP from transcript id
                info = fetch_transcript(self.transcript_id)
            if self.transcript_id is None:
                if self.gene_id is not None:
                    # REST LOOKUP from gene id
                    vals = fetch_transcript_info_from_gene_id(self.gene_id)
                    self.transcript_id, info = vals
                else:
                    raise LookupFailedException("No match for {}".
                                                format(self.gene_id))
        except LookupFailedException:
            # REST LOOKUP from symbol
            self.gene_id = fetch_gene_id_from_symbol(hugo_symbol)
            # REST LOOKUP from gene id
            vals = fetch_transcript_info_from_gene_id(self.gene_id)
            self.transcript_id, info = vals

        if lookup_seq:
            self._cds_seq = self.cds_seq

        self.chrom = info['seq_region_name']
        self.strand = info['strand']
        self.n_exons = len(info['Exon'])
        exon_starts = [i['start'] for i in info['Exon']]
        exon_ends = [i['end'] for i in info['Exon']]
        self.cdna_len = np.sum(np.array(exon_ends) - np.array(exon_starts) + 1)

        u5 = [i for i in info['UTR'] if i['object_type'] == 'five_prime_UTR']
        u3 = [i for i in info['UTR'] if i['object_type'] == 'three_prime_UTR']

        u5_cutoff, u3_cutoff = self._get_utr_cutoffs(u5, u3)

        self.cds_intervals = self._get_cds_intervals(
            self.strand, exon_starts, exon_ends, u5_cutoff, u3_cutoff)
        self.n_cds_intervals = len(self.cds_intervals)
        a, b = zip(*self.cds_intervals)
        self.cds_len = sum(np.array(b) - np.array(a) + 1)

        if lookup_seq:
            if self.cds_len != len(self.cds_seq):
                raise LengthMismatchException(
                    "CDS intervals don't match seq "
                    "length for {}".format(hugo_symbol))

        self.n_codons = int(self.cds_len / 3)

    def __repr__(self):
        return "<{classname} for {hugo}. {n_aa}aa, {cdslen}bp>".format(
            classname=self.__class__.__name__, hugo=self.hugo_symbol,
            n_aa=self.n_codons - 1, cdslen=self.cds_len)

    def _get_utr_cutoffs(self, u5, u3):
        u5_cutoff, u3_cutoff = None, None
        u5_starts = [i['start'] for i in u5]
        u5_ends = [i['end'] for i in u5]
        u3_starts = [i['start'] for i in u3]
        u3_ends = [i['end'] for i in u3]
        if self.strand == -1:
            if u5:
                u5_cutoff = min(u5_starts)
            if u3:
                u3_cutoff = max(u3_ends)
        else:  # + strand
            if u5:
                u5_cutoff = max(u5_ends)
            if u3:
                u3_cutoff = min(u3_starts)
        return u5_cutoff, u3_cutoff

    @property
    def cds_seq(self):
        """Lookup transcript, exon and UTR info from Ensembl Rest API."""
        if not self._cds_seq:
            seq_info = client.perform_rest_action(
                '/sequence/id/{0}'.format(self.transcript_id),
                params={'type': 'cds'})
            seq_str = seq_info['seq']
            if seq_str.startswith('N') or seq_str.endswith('N'):
                seq_str = seq_str.strip('N')
            seq = Seq(seq_str, IUPAC.unambiguous_dna)
            self._cds_seq = seq
            self._aa_seq = seq.translate()
        else:
            seq = self._cds_seq
        return seq

    @property
    def aa_seq(self):
        """Lookup transcript, exon and UTR info from Ensembl Rest API."""
        if not self._aa_seq:
            _ = self.cds_seq  # load cds_seq and aa_seq from ensembl api
        return self._aa_seq

    @staticmethod
    def _get_cds_intervals(strand, exon_starts, exon_ends,
                           u5_cutoff, u3_cutoff):
        """
        Get exon intervals after stripping out UTRs.
        for either strand, exons always listed in 5-3 order.
        start position always less than end position.
        """
        exon_intervals = list(zip(exon_starts, exon_ends))
        # force increasing coordinate exon order
        if strand == -1:
            exon_intervals.reverse()
            low_cutoff = u3_cutoff
            high_cutoff = u5_cutoff
        else:
            low_cutoff = u5_cutoff
            high_cutoff = u3_cutoff
        # remove low_cutoff section
        if low_cutoff is None:
            keep_intervals1 = exon_intervals[:]
        else:
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
        if high_cutoff is None:
            keep_intervals2 = keep_intervals1
        else:
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


def fetch_ids_for_hugo(hugo_symbol):
    """Get canonical transcript id and gene_id from ensembl."""
    transcript_id, gene_id = None, None
    ids_df = ensembl_df.loc[hugo_symbol]
    try:
        gene_id, transcript_id = ids_df
    except:
        pass
    return transcript_id, gene_id


def fetch_gene_id_from_symbol(symbol):
    """Lookup ensembl gene symbol from Ensembl REST API."""
    symbol_lookup = client.perform_rest_action(
        '/xrefs/symbol/human/{0}'.format(symbol),
        params={'external_db': 'HGNC', 'object_type': 'gene'})
    if not symbol_lookup:
        raise LookupFailedException('Failed to find match for {}'.
                                    format(symbol))
    # THERE IS AT LEAST ONE SYMBOL MATCH
    ids = [i['id'] for i in symbol_lookup]
    if len(symbol_lookup) == 1:
        gene_id = ids[0]
    else:
        data = {'ids': ids}
        d = client.perform_rest_action('/lookup/id', data_dict=data)
        good_ids = [i for i in d if d[i]['display_name'] == symbol]
        if not good_ids:
            raise LookupFailedException("No exact matches for {}"
                                        .format(symbol))
        else:  # use first matching id
            gene_id = good_ids[0]
    return gene_id


def fetch_transcript(transcript_id):
    """Lookup transcript, exon and UTR info from Ensembl Rest API."""

    transcript_info = client.perform_rest_action(
        '/lookup/id/{0}'.format(transcript_id),
        params={'object_type': 'transcript',
                'expand': 1, 'utr': 1})
    return transcript_info


def fetch_transcript_info_from_gene_id(gene_id):
    """ Look up transcript, exon and UTR info from Ensembl Rest API."""
    gene = client.perform_rest_action(
        '/lookup/id/{0}'.format(gene_id),
        params={'expand': '1', 'utr': 1})

    transcripts = [i for i in gene['Transcript'] if i['is_canonical'] == 1]
    if len(transcripts) != 1:
        raise LookupFailedException("Lookup for gene {} failed".format(
            gene_id))
    transcript_info = transcripts[0]
    transcript_id = transcript_info['id']
    return transcript_id, transcript_info


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
