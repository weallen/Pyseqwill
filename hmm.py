from scipy import stats
import numpy
import common
import table
import math

import pyximport; pyximport.install()
import cHMM

# K = number of hidden states
# M = number of input marks
class HMM:
    def __init__(self, M, K):
        self.M = M
        self.K = K
        # prob that will transition from state i to state j
        self.trans_probs = numpy.zeros((self.K, self.K))
        # prob that in state k has mark m
        self.emission_probs = numpy.zeros((self.K, self.M)) 
        self.start_probs = numpy.zeros(self.K)

    def random_init(self):
        self.trans_probs = numpy.random.random((self.K, self.K))
        self.start_probs = numpy.random.random(self.K) 
        self.emission_probs = numpy.random.random((self.K, self.M))

    def prob(self, obs):
        pass

    # Baum-Welch algorithm for determining trans and emit probs
    def train(self, obs):
        pass

    # AFter "Numerically Stable Hidden markov Model Implementation" by Tobias Mann
    def forward(self, obs):
        T = obs.shape[0]
        alpha = numpy.zeros((T, self.K))
        scale = numpy.zeros(T)
        alpha[0] = self.start_probs * self.all_obs_prob(obs[0])
        scale[0] = sum(alpha[0])
        alpha[0] /= scale[0]
#        alpha[0, :] = numpy.ones(self.K)
        for t in xrange(1, T): 
            obs_probs = self.all_obs_prob(obs[t])
            alpha[t] = numpy.array([sum(alpha[t-1] * self.trans_probs[:,j]) for j in range(0, self.K)]) * obs_probs
            s = sum(alpha[t, :])
            scale[t] = s
            alpha[t] = alpha[t, :]/s
        print alpha[T-1] 
        return (alpha, scale)
 
    def backward(self, obs, scale):
        T = obs.shape[0]
        beta = numpy.zeros((T, self.K))
        obs_probs = self.all_obs_prob(obs[T-1])
        beta[T - 1] = numpy.array([sum(self.trans_probs[j,:] * obs_probs) for j in range(0, self.K)]) / scale[T - 1]
        #beta[T - 1, :] = self.start_probs
        for t in xrange(T-2, -1, -1):
            obs_probs = self.all_obs_prob(obs[t+1])
            beta[t] = numpy.array([sum(beta[t+1] * obs_probs * self.trans_probs[j,:]) for j in range(0, self.K)]) / scale[t]
        #return sum(self.starts_probs * self.all_obs_prob(obs[0]) * beta[0])
        print beta[0]
        return beta
 
    def _count_transitions(self, obs):
        A = numpy.zeros((self.K, self.K))
        pass

    def viterbi_decode(self, obs):
        path = numpy.zeros(len(obs))
        pass
  
    # Returns an array of the likelihood of emitting obs for all hidden states  
    def all_obs_prob(self, obs):
        return numpy.array([self.obs_prob(s, obs) for s in range(0, self.K)])

    # Determines probaility distribution over set of observations 
    def obs_prob(self, state, obs):
        prob = 1.0
        for i in range(0,len(obs)):
            call = obs[i] # must be either 0 or 1
            emit_p = self.emission_probs[state, i] # prob of mark i in state
            if call == 1.0:
                prob *= emit_p
            else:
                prob *= 1 - emit_p
        return prob

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
    hmm = HMM(len(common.DATA_SETS), K) 
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
        alpha, scale = hmm.forward(obs)
        print "got to backward"
        beta = hmm.backward(alpha, scale)
    return hmm 
