from scipy import stats
import numpy
import common
import table

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
        self.trans_probs = numpy.random((self.K, self.K))
   
    def prob(self, obs):
        pass

    # Baum-Welch algorithm for determining trans and emit probs
    def train(self, obs):
        pass

    def forwards(self, obs):
        pass
 
    def backwards(self, obs):
        pass
    
    def _count_transitions(self, obs):
        A = numpy.zeros((self.K, self.K))
        pass

    def viterbi_decode(self, obs):
        path = numpy.zeros(len(obs))
        pass

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
 
def call_row(row, thresholds):
    calls = numpy.zeros(len(common.DATA_SETS))
    i = 0
    for name in common.DATA_SETS:
        if row[name] >= thresholds[name]:
            calls[i] = 1
        i += 1
    return calls

def chr_obs_seq(chr_cov, thresholds):
    obs = []
    for row in chr_cov.coverage:
        obs.append(call_row(row, thresholds))
    return obs

def find_empirical_means(tab):
    means = {}
    num_windows = len(tab.cols.start)
    for name in common.DATA_SETS:
        means[name] = sum(row[name] for row in tab)/num_windows
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
    means = find_empirical_means(tab)
    print means
    thresholds = find_threshold_vals(means)    
    print thresholds
    for chr in common.CHROMOSOMES:
        chr_cov = table.Chromosome(tab, chr)
        obs = chr_obs_seq(chr_cov, thresholds) 
   #     hmm.train(obs)
    return hmm


 
