#!/usr/bin/env python
import numpy as np
from sqlalchemy import text
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from ensembl_client import EnsemblRestClient, LookupFailedException
from . import NoIntervalsException, engine


class LengthMismatchException(Exception):
    pass


class CanonicalInfo(object):
    """Builds on ensembl api transcript info object.

    Provide hugo_symbol and (optionally) ensembl transcript and gene to cut down
    on lookups. Specify lookup_seq=True for cds and aa sequence lookups.
    """
    def __init__(self, hugo_symbol=None, transcript_id=None, ensembl_gene=None,
                 lookup_seq=False, server='http://grch37.rest.ensembl.org'):
        self.client = EnsemblRestClient(server=server)
        self.hugo_symbol = hugo_symbol
        self.transcript_id = transcript_id
        self.gene_id = ensembl_gene
        info = None

        try:
            # fetch_ids_for_hugo OR fetch_transcript could fail
            if transcript_id is None:
                transcript_gene = self.fetch_ids_for_hugo(hugo_symbol)
                self.transcript_id, self.gene_id = transcript_gene
            if self.transcript_id is not None:
                # REST LOOKUP from transcript id
                info = self.fetch_transcript(self.transcript_id)
            if self.transcript_id is None:
                if self.gene_id is not None:
                    # REST LOOKUP from gene id
                    vals = self.fetch_transcript_info_from_gene_id(self.gene_id)
                    self.transcript_id, info = vals
                else:
                    raise LookupFailedException("No match for {}".
                                                format(self.gene_id))
        except LookupFailedException:
            # REST LOOKUP from symbol
            self.gene_id = self.fetch_gene_id_from_symbol(hugo_symbol)
            # REST LOOKUP from gene id
            vals = self.fetch_transcript_info_from_gene_id(self.gene_id)
            self.transcript_id, info = vals
        
        if lookup_seq:
            self.cds_seq = self.get_cds_seq(self.transcript_id)
            self.aa_seq = self.cds_seq.translate()

        self.chrom = info['seq_region_name']
        self.strand = info['strand']
        self.n_exons = len(info['Exon'])
        exon_starts = [i['start'] for i in info['Exon']]
        exon_ends = [i['end'] for i in info['Exon']]
        self.cdna_len = np.sum(np.array(exon_ends) - np.array(exon_starts) + 1)

        u5 = [i for i in info['UTR'] if i['object_type'] == 'five_prime_UTR']
        u3 = [i for i in info['UTR'] if i['object_type'] == 'three_prime_UTR']

        u5_cutoff, u3_cutoff = self._get_utr_cutoffs(u5, u3)

        cds_intervals = self._get_cds_intervals(self.strand, exon_starts,
                                                exon_ends, u5_cutoff, u3_cutoff)
        self.cds_intervals = cds_intervals
        self.n_cds_intervals = len(cds_intervals)
        a, b = zip(*cds_intervals)
        self.cds_len = sum(np.array(b) - np.array(a) + 1)
        
        if lookup_seq:
            if self.cds_len != len(self.cds_seq):
                raise LengthMismatchException("CDS intervals don't match seq "
                    "length for {}".format(hugo_symbol))

        self.n_codons = self.cds_len/3

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

    def fetch_transcript(self, transcript_id):
        """Lookup transcript, exon and UTR info from Ensembl Rest API."""

        transcript_info = self.client.perform_rest_action(
            '/lookup/id/{0}'.format(transcript_id),
            params={'object_type': 'transcript',
                    'expand': 1, 'utr': 1})
        return transcript_info

    def fetch_gene_id_from_symbol(self, symbol):
        """Lookup ensembl gene symbol from Ensembl REST API."""
        symbol_lookup = self.client.perform_rest_action(
            '/xrefs/symbol/human/{0}'.format(symbol),
            params={'external_db': 'HGNC', 'object_type':'gene'})
        if not symbol_lookup:
            raise LookupFailedException('Failed to find match for {}'.\
                                        format(symbol))
        # THERE IS AT LEAST ONE SYMBOL MATCH
        ids = [i['id'] for i in symbol_lookup]
        if len(symbol_lookup) == 1:
            gene_id = ids[0]
        else:
            data = {'ids': ids}
            d = self.client.perform_rest_action('/lookup/id', data_dict=data)
            good_ids = [i for i in d if d[i]['display_name'] == symbol]
            if not good_ids:
                raise LookupFailedException("No exact matches for {}"
                                            .format(symbol))
            else:  # use first matching id
                gene_id = good_ids[0]
        return gene_id

    def fetch_transcript_info_from_gene_id(self, gene_id):
        """ Look up transcript, exon and UTR info from Ensembl Rest API."""
        gene = self.client.perform_rest_action(
            '/lookup/id/{0}'.format(gene_id),
            params={'expand': '1', 'utr': 1})

        transcripts = [i for i in gene['Transcript'] if i['is_canonical'] == 1]
        if len(transcripts) != 1:
            raise LookupFailedException("Lookup for gene {} failed".format(
                gene_id))
        transcript_info = transcripts[0]
        transcript_id = transcript_info['id']
        return transcript_id, transcript_info

    def get_cds_seq(self, transcript_id):
        """Lookup transcript, exon and UTR info from Ensembl Rest API."""
        seq_info = self.client.perform_rest_action(
            '/sequence/id/{0}'.format(transcript_id),
            params={'type': 'cds'})
        seq_str = seq_info['seq']
        if seq_str.startswith('N') or seq_str.endswith('N'):
            seq_str = seq_str.strip('N')
        seq = Seq(seq_str, IUPAC.unambiguous_dna)
        return seq

    @staticmethod
    def fetch_ids_for_hugo(hugo_symbol):
        """Get canonical transcript id and gene_id from ensembl."""
        con, rows, transcript_id, gene_id = None, None, None, None
        cmd1 = "SELECT `canonical_transcript`, `ensembl_gene_id`,  hg19_valid "\
               "FROM refs.`ensembl_canonical` WHERE hugo_symbol = :hugo;".\
            format(hugo_symbol)
        cur = engine.execute(text(cmd1), hugo=hugo_symbol)
        if cur.rowcount != 1:
            raise LookupFailedException("Non-single row for transcript {}"
                                        .format(hugo_symbol))
        rows = cur.fetchall()
        # rows is [[transcript,gene],[transcript,gene],...]
        for row in rows:
            hg19_valid = row[2]
            if hg19_valid == 1:
                transcript_id = str(row[0])
            else:
                transcript_id = None
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