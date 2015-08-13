#!/usr/bin/env python
import MySQLdb as mdb
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from numpy import cumsum, array, concatenate, nonzero
import subprocess
import os
import pdb
from lookup_hg19 import test_cpg, lookup_hg19

dbvars = {'host': 'localhost', 'db': 'CancerDB',
          'read_default_file': "~/.my.cnf"}

bed_bin_dir = os.path.join(os.path.expanduser('~'),'bin','bed')  # bedtools

class NoIntervalsException(Exception):
    pass


class LookupFailedException(Exception):
    pass


class GeneSeq():
    """Holds exon CDS info for obtaining transcript coordinates."""

    def __init__(self, refseq_NM):
        """Takes string inputs."""
        self.refseq = refseq_NM
        self.exonSet = ExonSet(self.refseq)
        # self.seq = self._fetchMrnaSeq()
        self.cds_seq = self._fetch_cds_seq(self.exonSet)
        self.hugo = self.exonSet.hugo
        utr_5p_len = self.exonSet.utr_5p_len
        utr_3p_len = self.exonSet.utr_3p_len
        cds_len = self.exonSet.cds_len
        # if len(self.seq) != (utr_5p_len + utr_3p_len + cds_len):
        #     raise Exception(
        #         "Seq length {} doesn't match utr and cds lengths {} for {}.".format(
        #             refseq_NM))
        # cds_seq = self.seq[utr_5p_len:-utr_3p_len]
        # cds_seq = Seq(cds_seq, IUPAC.unambiguous_dna)
        # self.cds_seq = cds_seq
        self.aa_seq = self.cds_seq.translate()
        self.n_codons = cds_len / 3

    def __repr__(self):
        return "<{classname} for {hugo} ({refseq}). n_aa={n_aa}, cds_len={cdslen}>".format(
            classname=self.__class__.__name__, hugo=self.hugo,
            refseq=self.refseq,
            n_aa=self.n_codons - 1, cdslen=self.exonSet.cds_len)

    def _fetchMrnaSeq(self, hgVersion='hg19'):
        """NOT USED."""
        rows = None
        seq = None
        # GET pway_size
        cmd1 = """SELECT seq FROM {}.knownGeneTxMrna WHERE refseq_NM = {!r};""".format(
            hgVersion, self.refseq)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd1)
            rowCount = cur.rowcount
            if rowCount != 1:
                raise LookupFailedException(
                    "Non-single row found for transcript {}".format(
                        self.refseq))
            rows = cur.fetchall()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        # rows is [[id,name],[id,name],...]
        for row in rows:
            seq = str(row[0])
        return seq

    def _fetch_cds_seq(self, exonSet):
        seq_list = list()
        chrom = exonSet.chrom.lstrip('chr')
        for interval in exonSet.cds:
            seq_list.append(lookup_hg19(chrom, *interval))
        if exonSet.strand == '-':
            seq_list = seq_list[::-1]
            seq_list = [i[::-1] for i in seq_list]
        seq = ''.join(seq_list)
        seq = Seq(seq, IUPAC.unambiguous_dna)
        if exonSet.strand == '-':
            seq = seq.complement()
        return seq


class ExonSet():
    """Holds exon CDS info for obtaining transcript coordinates."""

    def __init__(self, refseq_NM):
        """Takes string inputs."""
        self.refseq = refseq_NM
        (strand, cdsStart, cdsEnd, exonStarts, exonEnds, chrom,
         hugo) = self._fetchTranscriptInfo()
        cdsStart += 1  # add 1 to get 1-based coords
        self.strand = strand
        self.cds_start = cdsStart
        self.cds_end = cdsEnd
        self.chrom = chrom
        self.hugo = hugo
        self.cds = None
        self.cds_len = None
        self.utr_5p = None
        self.utr_3p = None
        self.utr_5p_len = None
        self.utr_3p_len = None

        # exonStarts = array([int(i) for i in exonStarts.rstrip(',').split(',')]) + 1  # add 1 to get 1-based coords
        # exonEnds = array([int(i) for i in exonEnds.rstrip(',').split(',')])
        exonStarts = [int(i) + 1 for i in exonStarts.rstrip(',').split(
            ',')]  # add 1 to get 1-based coords
        exonEnds = [int(i) for i in exonEnds.rstrip(',').split(',')]
        if len(exonStarts) != len(exonEnds):
            raise Exception("Exon starts and ends arrays are not same size.")
        (utr_a, self.cds, utr_b) = self.divide_exons(exonStarts, exonEnds)
        # assign utr_5p and utr_3p according to strand.
        if strand == '+':
            self.utr_5p = utr_a
            self.utr_3p = utr_b
        elif strand == '-':
            self.utr_5p = utr_b
            self.utr_3p = utr_a
        else:
            raise Exception(
                "Invalid strand provided for {}".format(self.refseq))
        self.n_exons_cds = len(self.cds)
        self.cds_len = sum([i[1] - i[0] + 1 for i in self.cds])
        self.utr_5p_len = sum([i[1] - i[0] + 1 for i in self.utr_5p])
        self.utr_3p_len = sum([i[1] - i[0] + 1 for i in self.utr_3p])

        if self.cds_len % 3 != 0:
            print("{} has a coding length of {} - {} mod 3.".format(self.refseq,
                                                                    self.cds_len,
                                                                    self.cds_len % 3))
            # self.intervals = zip(exonStarts, exonEnds)
            # self.interval_lengths = tuple(exonEnds - exonStarts + 1)

    def divide_exons(self, exonStarts, exonEnds):
        """Assign (start,end) coordinate tuples to lists for utr_a, cds,
        and utr_b."""
        cdsStart = self.cds_start
        cdsEnd = self.cds_end
        # Three lists will hold tuples of start,end coordinates. utr_a is utr with lowest coords on genome ref.
        utr_a = list()
        cds = list()
        utr_b = list()
        # CDS overlap tests. iterate over exons, assigning to utrs or cds.
        for ind in xrange(len(exonStarts)):
            exStart = exonStarts[ind]
            exEnd = exonEnds[ind]
            # if interval before cds
            if exEnd < cdsStart:
                utr_a.append((exStart, exEnd))
            # if interval overlaps cdsStart
            elif exStart < cdsStart and exEnd >= cdsStart:
                utr_a.append((exStart, cdsStart - 1))
                cds.append((cdsStart, exEnd))
            # if interval contained within cds
            elif exStart >= cdsStart and exEnd <= cdsEnd:
                cds.append((exStart, exEnd))
            # if interval overlaps cdsEnd
            elif exStart <= cdsEnd and exEnd > cdsEnd:
                cds.append((exStart, cdsEnd))
                utr_b.append((cdsEnd + 1, exEnd))
            # interval comes after cds
            elif exStart > cdsEnd:
                utr_b.append((exStart, exEnd))
            else:
                raise Exception("Mysteriously failed cds overlap tests")
        return (utr_a, cds, utr_b)

    def _get_CDS_coord_plus(self, position):
        """Get coordinate, indexing from 0, in CDS, positive strand."""
        position = int(position)  # ensure integer
        # get interval index containing position
        use_interval = None
        found_interval = False
        for ind in xrange(self.n_exons_cds):
            interval = self.cds[ind]
            if position >= interval[0] and position <= interval[1]:
                found_interval = True
                use_interval = ind
                interval_pos = position - interval[0]
                break
        if not found_interval:
            raise NoIntervalsException(
                "No interval found for position {} in {}".format(position,
                                                                 self.refseq))
        cds_lengths = [i[1] - i[0] for i in self.cds]
        prev_length = sum(cds_lengths[:use_interval])
        coord = interval_pos + prev_length  # e.g. prevlen of 5 and interval_pos of 0 should give coord=5
        return coord

    def get_CDS_coord(self, position):
        """Get coordinate, indexing from 0, in CDS, on gene-specific strand."""
        temp_coord = self._get_CDS_coord_plus(position)
        if self.strand == '+':
            coord = temp_coord
        elif self.strand == '-':
            coord = self.cds_len - 1 - temp_coord
        return coord

    def _get_hg_coord_plus(self, cds_ind):
        """Get human genome coordinate corresponding to cds_coordinate (indexed from 0)."""
        cds_ind = int(cds_ind)  # ensure integer

        cds_lens = [i[1] - i[0] + 1 for i in self.cds]
        cumul = cumsum(array(cds_lens))
        cumul_zero = concatenate((array([0]), cumul), axis=1)
        # position is in first interval where cumsum of previous intervals is GREATER than cds_index
        interval_ind = nonzero(cumul > cds_ind)[0][0]
        subinterval_ind = cds_ind - cumul_zero[interval_ind]
        hg_pos = self.cds[interval_ind][0] + subinterval_ind
        return hg_pos

    def get_hg_coord(self, cds_ind):
        """Get human genome coordinate corresponding to cds_coordinate on transcript-specific strand."""
        if self.strand == '+':
            hg_pos = self._get_hg_coord_plus(cds_ind)
        elif self.strand == '-':
            cds_ind_plus = self.cds_len - 1 - cds_ind  # get cds index on + strand
            hg_pos = self._get_hg_coord_plus(cds_ind_plus)
        return hg_pos

    def _fetchTranscriptInfo(self, hgVersion='hg19'):
        """Get pathway ids containing genes in (possibly empty) interest set."""
        rows = None
        strand = None
        cdsStart = None
        cdsEnd = None
        exonStarts = None
        exonEnds = None
        # GET pway_size
        cmd1 = """SELECT strand, cdsStart, cdsEnd, exonStarts, exonEnds, chrom, hugo_symbol
         FROM {}.refFlat WHERE `name` = {!r};""".format(hgVersion, self.refseq)
        try:
            con = mdb.connect(**dbvars)
            cur = con.cursor()
            cur.execute(cmd1)
            rowCount = cur.rowcount
            if rowCount != 1:
                raise Exception("Non-single row found for transcript {}".format(
                    self.refseq))
            rows = cur.fetchall()
        except mdb.Error as e:
            print "Error %d: %s" % (e.args[0], e.args[1])
        finally:
            if con:
                con.close()
        # rows is [[id,name],[id,name],...]
        for row in rows:
            strand = str(row[0])
            cdsStart = int(row[1])
            cdsEnd = int(row[2])
            exonStarts = str(row[3])
            exonEnds = str(row[4])
            chrom = str(row[5])
            hugo = str(row[6])
        return (strand, cdsStart, cdsEnd, exonStarts, exonEnds, chrom, hugo)


class TranscriptSiteInfo():
    """Populated by append_category_beds."""

    def __init__(self, hugo, refseq, chrom):
        self.hugo = hugo
        self.refseq = refseq
        self.chrom = chrom
        self.site_list = list()

    # def add_site(self, pos, is_at, is_cg, is_cpg, has_silent, has_nonsilent):
    # """Add site category info to collection that will include whole transcript."""
    # self.site_list.append((pos, is_at, is_cg, is_cpg, has_silent, has_nonsilent))
    def add_site(self, pos, has_silent, has_nonsilent):
        """Add site category info to collection that will include whole
        transcript."""
        self.site_list.append((pos, has_silent, has_nonsilent))

    # def write_category_bed(self, file_obj, is_at, is_cg, is_cpg, has_silent, has_nonsilent):
    def write_category_bed(self, file_obj, write_silent, write_nonsilent):
        """Write to file sites that match desired category booleans."""
        for site in self.site_list:
            # if (site[1] is is_at and site[2] is is_cg and site[3] is is_cpg and
            #             site[4] is has_silent and site[5] is has_nonsilent):
            if (write_silent and site[1]) or (write_nonsilent and site[2]):
                pos_str = "{}\t{}\t{}\t{}\n".format(self.chrom, site[0] - 1,
                                                    site[0], self.hugo)
                file_obj.write(pos_str)


def append_category_beds(refseq_NM, silent_path, nonsilent_path):
    try:
        g = GeneSeq(refseq_NM)
    except LookupFailedException as e:
        print e.message
        return
    # print(g.cds_seq)
    alt_dict = dict(A='CGT', C='AGT', G='ACT', T='ACG')
    site_info = TranscriptSiteInfo(g.hugo, refseq_NM, g.exonSet.chrom)
    # bool_dict = {True:'1', False:'0'}
    # with open(g.hugo + "_silent.bed",'w') as out_silent:
    # with open(g.hugo + "_nonsilent.bed",'w') as out_nonsilent:
    # out_silent.write('chr\tstart\tend\n')
    # out_nonsilent.write('chr\tstart\tend\n')
    for codon_ind in xrange(g.n_codons):
        # identify previous and next base (prv_base, nxt_base)
        bases3 = g.cds_seq[codon_ind * 3:codon_ind * 3 + 3]
        aa_orig = g.aa_seq[codon_ind]
        # for each of 3 bases
        for subind in xrange(3):
            has_silent = False
            has_nonsilent = False
            cds_ind = codon_ind * 3 + subind
            hg_pos = g.exonSet.get_hg_coord(cds_ind)
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
            for alt in alt_dict[bases3[subind]]:
                newbases = bases3.tomutable()
                newbases[subind] = alt
                new_aa = newbases.toseq().translate()
                if new_aa[0] == aa_orig:
                    has_silent = True
                else:
                    has_nonsilent = True
                    # if has_nonsilent and has_silent:
                    # 	break  # out of alternative iteration
            site_info.add_site(hg_pos, has_silent, has_nonsilent)
    with open(silent_path, 'a') as out:
        site_info.write_category_bed(out, True, False)
    with open(nonsilent_path, 'a') as out:
        site_info.write_category_bed(out, False, True)
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


def condense_beds_silent_nonsilent(silent_path, nonsilent_path):
    """Sort and merge both silent and nonsilent raw bed files for gene."""
    for raw_path in [silent_path, nonsilent_path]:
        sorted_path = raw_path + "_sorted.bed"
        condense_path = raw_path + "_final.bed"

        cmd1 = ["sort", "-k1,1", "-k2,2n", raw_path, "-o", sorted_path]
        cmd2 = [os.path.join(bed_bin_dir,"mergeBed"), "-i", sorted_path, "-c", "4", "-o", "distinct"]
        subprocess.check_call(cmd1)
        # with open(condense_path, 'w') as file:
        #     subprocess.check_call(cmd2, stdout=file)
        with open(os.devnull, "r") as fnullin:
            with open(condense_path, "w") as outfile:
                p = subprocess.Popen(cmd2, stdin=fnullin, stdout=outfile)
                p.wait()
        os.remove(raw_path)
        os.remove(sorted_path)


def get_beds_for_refseq_list(refseq_list, silent_path, nonsilent_path):
    """Takes list of refseq_NM strings and creates silent and nonsilent bed
    files."""
    for path in [silent_path, nonsilent_path]:
        if os.path.exists(path):
            os.remove(path)

    for refseq in refseq_list:
        append_category_beds(refseq, silent_path, nonsilent_path)
    condense_beds_silent_nonsilent(silent_path, nonsilent_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument("refseq", help="Refseq_NM.")
    parser.add_argument("refseq_file", help="File with list of Refseq_NM ids.",
                        type=file)
    args = parser.parse_args()
    # refseq = args.refseq
    # append_category_beds(refseq)
    refseq_file = args.refseq_file
    refseq_list = list()
    try:
        for line in refseq_file:
            refseq_list.append(line.strip().split('\t')[1])
    finally:
        refseq_file.close()

    silent_path = "/Users/Stephen/Downloads/roi_silent"
    nonsilent_path = "/Users/Stephen/Downloads/roi_nonsilent"
    get_beds_for_refseq_list(refseq_list, silent_path, nonsilent_path)

    # g = GeneSeq(refseq)
    # print g












