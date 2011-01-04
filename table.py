import common
import tables
import data
import numpy as np

class GenomeWindow(tables.IsDescription):
    idx = tables.Int32Col() # absolute index of window
    chr = tables.Int8Col() # chromosome 
    start = tables.Int32Col() # end of window
    end = tables.Int32Col() # start of window
    chr_pos = tables.Int32Col()

class Chromosome:
    def __init__(self, tab, chr):
        self.tab = tab
        self.chr = common.CHR_TO_NUM[chr]
        self.coverage = [r for r in self.tab.readWhere("chr == %d" % self.chr)]
        sorted(self.coverage, key=lambda x: x["start"])
    

class SeqData:
    def __init__(self, fname):
        self.h5file = tables.openFile(fname, mode="r")
        self.windows = self.h5file.root.windows
        self.chr_bounds = {}
        for chr in common.CHROMOSOMES:
            self.chr_bounds[chr] = self._find_chr_bounds_idx(chr)
        self.data = {}
        for d in common.DATA_SETS:
            self.data[d] = self.h5file.root.seqdata._f_getChild(d).read()

    # returns bounds [start, end)
    def _find_chr_bounds_idx(self, chr):
        if chr in common.CHROMOSOMES:
            chr_idx = common.CHR_TO_NUM[chr]
            chr_start = 0
            chr_stop = 0
            for row in self.windows.where("chr == %i" % (chr_idx)):
                idx = row['idx']
                if row['chr_pos'] == 0: 
                    chr_start = idx
                if idx > chr_stop:
                    chr_stop = idx
            return (chr_start, chr_stop+1)
        return None

    def get_data_by_chr(self, dataset, chr):
        chr_start, chr_end = self.chr_bounds[chr]
        return self.data[dataset][chr_start:chr_end]
    
    def get_all_data_by_chr(self, chr):
        chr_len = self.chr_bounds[chr][1] - self.chr_bounds[chr][0]
        all = np.zeros((chr_len, len(common.DATA_SETS)))
        i = 0
        for d in common.DATA_SETS: 
            all[:, i] = self.get_data_by_chr(d, chr)
            i += 1
        return all
def load_seq_data():
    return SeqData(common.DATA_PATH+"all_data.h5")

def import_seq_data():
    genome = data.Genome("/gpfs/runtime/bioinfo/bowtie/indexes/mm9_all_folded_with_chrID.fa")
    print "Loaded genome"
    h5file = tables.openFile(common.DATA_PATH+"all_data.h5", mode="w", title="hme and chip data")
    create_window_table(h5file, genome)
    load_data_into_seq_arrays(h5file, genome)
    h5file.close()

def create_window_table(h5file, genome):
    print "Creating window table"
    tab = h5file.createTable(h5file.root, 'windows', GenomeWindow, "Sliding windows across genome")
    parts = genome.partition(common.WINDOW_SIZE)
#    num_parts = genome.num_parts(common.WINDOW_SIZE)
#    genome_array = np.zeros((num_parts, 5), dtype=np.int32) 
    i = 0
    window = tab.row
    for chr_idx in range(20):
        rngs = parts[chr_idx] 
        print "Chr", chr_idx+1
        n = 0
        for rng in rngs:
            window['idx'] = i
            window['chr' ] = chr_idx+1
            window['start'] = rng[0]
            window['end'] = rng[1]
            window['chr_pos'] = n
            window.append()
#            genome_array[i] = np.array([i, chr_idx, rng[0], rng[1], n]) 
            n += 1
            i += 1
    tab.flush()
    #h5file.createArray(h5file.root, 'windows', genome_array)


def load_data_into_seq_arrays(h5file, genome):
    group = h5file.createGroup(h5file.root, "seqdata", 'Seq data')
    for d in common.DATA_SETS:
        print "Loading data set ", d
        cov = data.Coverage(common.DATA_PATH + d + ".bam", genome)
        cov_arr = np.zeros(len(h5file.root.windows), dtype=np.int32)
        #genome_array = h5file.getNode("/", "windows").read()
        for window in h5file.root.windows:
            i = window['idx']
            chr_idx = window['chr']-1
            chr_pos = window['chr_pos']
            cov_arr[i] = cov.pileup[chr_idx][chr_pos]
        h5file.createArray(group, d, cov_arr)
