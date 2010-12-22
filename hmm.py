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

def chr_obs_seq(chr_cov, thresholds):
    obs = numpy.zeros((len(chr_cov.coverage), len(common.DATA_SETS)), int)
    for i in xrange(0,len(chr_cov.coverage)):
        obs[i] = call_row(chr_cov.coverage[i], thresholds)
    return obs

def find_empirical_means(tab):
    means = {}
    num_windows = float(sum(1 for row in tab))
    for name in common.DATA_SETS:
        means[name] = round(sum(row[name] for row in tab)/num_windows)
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

def train_hmm(tab, K):
    hmm = cHMM.HMM(len(common.DATA_SETS), K) 
    hmm.random_init()
    means = find_empirical_means(tab)
    print "MEANS ", means
    thresholds = find_threshold_vals(means)    
    print "THRESHOLDS ", thresholds
    A = numpy.zeros(hmm.M)
    E = numpy.zeros(hmm.M)
    for chr in common.CHROMOSOMES:
        print chr
        chr_cov = table.Chromosome(tab, chr)
        obs = chr_obs_seq(chr_cov, thresholds)
        print obs.shape[0]
        alpha, scale = cHMM.forward(hmm, obs)
        print "got to backward"
        beta = cHMM.backward(hmm, alpha, scale)
    return hmm 
