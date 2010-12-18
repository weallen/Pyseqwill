#!/usr/bin/env python

import common
import pysam
import tables
import features
import math
import itertools

class GenomePostion(tables.IsDescription):
    chr = StringCol(5)
    pos = UInt32Col()
    src = StringCol(16)
    type = StringCol(16)
    val = UInt32Col()
    
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
    f = pysam.Sam.file(bam_file_name, "rb")
    return f

class Genome:
    def __init__(self, fasta):
        self.fasta = fasta
        self.sizes = {}
        for chr in common.CHROMOSOMES:
            self.sizes = _chromosome_size(chr)
         
    def partition(self, part_size):
        chr_ranges = {}
        for chr in common.CHROMOSOMES:
            chr_ranges[chr] = []
            chr_size = self.sizes[chr]
            for i in xrange(0, chr_size, part_size):
                if i + part_size < chr_size:
                    chr_ranges[chr].append((i,i+part_size))
                else:
                    chr_ranges[chr].append((i, chr_size))
        return chr_ranges

    def _chromosome_size(self, chr):
        return len(self.fasta.fetch(chr))
    
    
class Coverage:
    def __init__(self, bam, genome):
        self.genome = genome
        self.part_size = 200
        self.ranges = genome.partition(self.part_size)
        self.bam = bam 
        for chr, parts in self.ranges.iteritems():
            num_intervals = len(parts)
            self.pileup[chr] = [i for i in itertools.repeat(0, num_intervals)]
         

    def compute_tag_overlap(self):
        pileup = {}
        for chr, rng in self.ranges.iteritems():
            p = self.bam.pileup(chr, rng[0], rng[1])
            for col in p:
                start = col.pos + 350
                chr_len = self.genome.sizes[chr]
                bin = int(math.floor(start / chr_len))
                pileup[chr][bin] += col.n
        return pileup

