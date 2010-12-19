#!/usr/bin/env python

import common
import pysam
import tables
import features
import math
import itertools
import numpy

       
class GFFFeature:
    def __init__(self, seqname, source, feature, start,\
                    end, score, strand, frame, attr):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attr = attr
 

def parse_gtf_attr(attr):
    attrs = attr.split(";")
    # the next line takes 'gene_id "blah"' -> blah
    gene_id = attrs[0].split(" ")[1][1:-1]
    transcript_id = attrs[1].split(" ")[1][1:-1]
    return {'gene_id': gene_id, 'transcript_id' : transcript_id}

def load_kg_gtf(gtf_file_name):
    f = pysam.TabixFile(gtf_file_name)
    gtf = f.fetch(parser=pysam.asGTF())
    feats = []
    for row in gtf:
        attr = parse_gtf_attr(row.attributes)
        currfeat = GFFFeature(row.seqname, row.source, row.feature,\
                                int(row.start), int(row.end), score, \
                                strand, frame, attr)
        feats.append(currfeat)
    return feats

def load_kgxref(kgxref_file_name):
    pass

def build_kg_transcripts(feats):
    tx = {}
    for feat in feats:
        txid = feat.attr["transcript_id"]
        if txid not in tx:
            tx[txid] = Transcript(txid)
        if feat.feature == "exon":
            tx[txid].exons.append(GRange(feat.seqname, feat.start, \
                                        feat.end, feat.strand))
            if not tx[txid].tx_range:
                tx[txid].tx_range = GRange(feat.seqname, feat.start, feat.end, feat.strand)
            else:
                if tx[txid].tx_range.start > feat.start:
                    tx[txd].tx_range.start = feat.start
                if tx[txid].tx_range.end < feat.end:
                    tx[txid].tx_range.end = feat.end
        elif feat.feature == "start_codon":
            if not tx[txid].cds_range:
                tx[txid].cds_range = GRange(feat.seqname, feat.start, 0, feat.strand)
            else:
                tx[txid].cds_range.start = feat.start
        elif feat.feature == "stop_codon":
            if not tx[txid].cds_range:
                tx[txid].cds_range = GRange(feat.seqname, 0, feat.end, feat.strand)
            else:
                tx[txid].cds_range.end = feat.end
    return tx

def build_kg_genes(trans, kgxref):
    genes = {}
    for sym, txid in kgxref:
        if not sym in genes:
            genes[sym] = Gene(sym)
        if txid in trans:
            genes[sym].add_transcript(trans[txid])
    return genes

def load_genome(genome_fasta):
    f = pysam.Fastafile(genome_fasta)
    return f

def load_bam_file(bam_file_name):
    f = pysam.Samfile(bam_file_name, "rb")
    return f

class Genome:
    def __init__(self, fasta):
        self.fasta = load_genome(fasta)
        self.sizes = {}
        for chr in common.CHROMOSOMES:
            self.sizes[chr] = self._chromosome_size(chr)
         
    def partition(self, part_size):
        parts = {}
        for chr in common.CHROMOSOMES:
            parts[chr] = PartitionIterator(self.sizes[chr], part_size) 
        return parts

    def _chromosome_size(self, chr):
        return len(self.fasta.fetch(chr))

class PartitionIterator:
    def __init__(self, chr_size, part_size):
        self.chr_size = chr_size
        self.part_size = part_size
        self.i = 0

    def next(self):
        start = self.i
        if start == self.chr_size:
            raise StopIteration
        if self.i + self.part_size < self.chr_size:
            self.i += self.part_size
        else:
            self.i = self.chr_size
        return (start, self.i)  
 
    def __iter__(self):
        return self

    def __len__(self):
        return int(math.floor(self.chr_size / self.part_size))
 
class Coverage:
    def __init__(self, bam, genome):
        self.genome = genome
        self.part_size = common.WINDOW_SIZE
        self.bam = load_bam_file(bam)
        self.pileup = {}
        for chr in common.CHROMOSOMES:
            num_intervals = int(math.floor(self.genome.sizes[chr] / self.part_size))
            self.pileup[chr] = numpy.zeros(num_intervals+1)
        self._compute_tag_overlaps()
        self.bam.close()

    # TODO Check with brad about which way to extend start
    def _compute_tag_overlaps(self):
        for chr in common.CHROMOSOMES:
            chr_len = self.genome.sizes[chr]
            reads = self.bam.fetch(chr, 0, chr_len)
            for read in reads:
                start = read.pos + 350
                if start > chr_len:
                    start = chr_len
                bin = int(math.floor(start / self.part_size))
                self.pileup[chr][bin] += 1
