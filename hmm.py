from scipy import stats
import numpy
import common
import table
import math

#import pyximport; pyximport.install()
import cHMM

# HMM Training stuff

def train_hmm(sd, K):
    hmm = cHMM.HMM(len(common.DATA_SETS), K) 
    hmm.random_init()
    A = numpy.zeros(hmm.M)
    E = numpy.zeros(hmm.M)
    for chr in common.CHROMOSOMES:
        print chr
        chr_cov = sd.get_all_data_by_chr(chr)
        print "Got chromosome coverage"
        print "Setting obs"
        hmm.set_obs(chr_cov)
        print "Got to trainng"
        hmm.train()
        print "got to backward"
    return hmm 
