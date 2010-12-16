
class GRange:
    def __init__(self, chr, start, end, strand):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand

    def __lt__(self, other):
        pass

    def __gt__(self, other):
        pass

    def __eq__(self, other):
        return (self.chr == other.chr) and (self.start == other.start)\
                and (self.end == other.end) and (self.strand == strand)

class Transcript:
    def __init__(self, txid):
        self.txid = txid
        self.tx_range = None
        self.cds_range = None
        self.exons = exons

class Gene:
    def __init__(self, sym):
        self.sym = sym
        self.range = None
        self.tx = []
        self.exons = []

    def add_transcript(tx):
        self.tx.append(tx)
        if not self.range:
            self.range = tx.tx_range
        else:
            if self.range.start < tx.tx_range.start:
                self.range.start = tx.tx_range.start
            if self.range.end > tx.tx_range.end:
                self.range.end = tx._tx_range.end
        uniq = True
        for txex in tx.exons:
            for myex in self.exons:
                if myex == txex:
                    uniq = False
                    break
            if uniq:
                self.exons.append(txex)
 
