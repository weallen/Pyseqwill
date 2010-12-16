#!/usr/bin/env python

import pysam

import features

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

def load_bam_file(bam_file_name):
    f = pysam.Sam.file(bam_file_name, "rb")
    return f
