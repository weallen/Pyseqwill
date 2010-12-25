from scipy import stats
import numpy
import common
import table
import math

#import pyximport; pyximport.install()
import cHMM

def is_over_thresh(val, thresh):
    if val >= thresh:
        return 1
    return 0

def call_row(row, thresholds):
    calls = numpy.array([is_over_thresh(row[name], thresholds[name]) for name in common.DATA_SETS])
    return calls

# chr_cov is 2D (num datasets x num windows) matrix of window counts
def chr_obs_seq(chr_cov, thresholds):
    obs = numpy.zeros((len(chr_cov), len(common.DATA_SETS)), dtype=numpy.int32)
    for i in xrange(0,chr_cov):
        obs[i] = call_row(chr_cov[i], thresholds)
    return obs

def find_empirical_means(sd):
    means = {}
    for name in common.DATA_SETS:
        means[name] = round(numpy.mean(sd.data[name]))
    return means

def find_threshold_vals(means):
    thresholds = {}
    for name in common.DATA_SETS:
        for i in range(1,15):
            mean = means[name]
            if stats.poisson.pmf(i, mean) < 1E-4: 
                thresholds[name] = i
                break
    return thresholds

# HMM Training stuff

def train_hmm(sd, K):
    hmm = cHMM.HMM(len(common.DATA_SETS), K) 
    hmm.random_init()
    means = find_empirical_means(sd)
    print "MEANS ", means
    thresholds = find_threshold_vals(means)    
    print "THRESHOLDS ", thresholds
    A = numpy.zeros(hmm.M)
    E = numpy.zeros(hmm.M)
    for chr in common.CHROMOSOMES:
        print chr
        chr_cov = sd.get_all_by_chr(chr)
        obs = chr_obs_seq(chr_cov, thresholds)
        print alpha[len(obs)-1]
        print "got to backward"
        beta = hmm.backward(scale)
        print beta[0]
    return hmm 
