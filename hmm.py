from scipy import stats
import numpy
import common
import table
import math

LOGZERO = float('nan')

def eln(x):
    if x == 0:
        return LOGZERO
    elif x > 0:
        return math.log(x)
    else:
        except ValueError

def eexp(x):
    if math.isnan(x):
        return 0
    return math.exp(x)

def elnsum(x, y):
    if math.isnan(x) or math.isnan(y):
        if math.isnan(x):
            return y
        return x
    if x > y:
        x + eln(1 + math.exp(y - x)
    return y + eln(1 + math.exp(x - y)
            
def elnproduct(x, y):
    if math.isnan(x) or math.isnan(y):
        return LOGZERO
    return x + y

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

    # AFter "Numerically Stable Hidden markov Model Implementation"
    def forward(self, obs):
        T = obs.shape[0]
        alpha = numpy.zeros((T, self.K)) 
#        alpha[0] = numpy.ones(self.K)
        alpha[0] = self.start_probs * self.all_obs_prob(obs[0])
        for t in xrange(1, T): 
            obs_probs = self.all_obs_prob(obs[t])
            s = sum(obs_probs) 
            for i in range(0, self.K):              
                alpha[t, i] = sum(self.trans_probs[:,i] * alpha[t-1]) 
            alpha[t] *= obs_probs
        if t % 1000 ==0 :
            print t, alpha[t]
        return sum(alpha[T-1] * self.start_probs)
 
    def backward(self, obs):
        T = obs.shape[0]
        beta = numpy.zeros((T, self.K))
        obs_probs = self.all_obs_prob(obs[T-1])
        beta[T - 1] = numpy.array([sum(self.trans_probs[j,:] * obs_probs) for j in range(0, self.K)]) 
        for t in xrange(T-1, 0, -1):
            obs_probs = self.all_obs_prob(obs[t])
            beta[t-1] = numpy.array([sum(beta[t] * obs_probs * self.trans_probs[j,:]) for j in range(0, self.K)])
        return sum(self.starts_probs * self.all_obs_prob(obs[0]) * beta[0])
    
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
    for chr in common.CHROMOSOMES:
        chr_cov = table.Chromosome(tab, chr)
        obs = chr_obs_seq(chr_cov, thresholds) 
        print hmm.forward(obs)
    return hmm 
