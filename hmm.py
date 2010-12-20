from scipy import stats
import numpy
import common

# K = number of hidden states
# M = number of input marks
class HMM:
    def __init__(self, obs, means, thresholds, K):
        self.obs = obs # table of windows
        self.means = means
        self.thresholds = thresholds
        self.M = len(common.DATA_SETS)
        self.K = K
        # prob that will transition from state i to state j
        self.trans_probs = numpy.zeros((self.K, self.K))
        # prob that in state k has mark m
        self.emission_probs = numpy.zeros((self.K, self.M)) 
        self.start_probs = numpy.zeros(self.K)

    def random_init(self):
        self.trans_probs = numpy.random((self.K, self.K))

    def obj_prob(self):
        pass

    def row_obs_prob(self, state, row):
        calls = self._call_row(row)
        prob = 1.0
        for i in range(0,len(calls)-1):
            call = calls[i] # must be either 0 or 1
            emit_p = self.emission_probs[state, i] # prob of mark i in state
            if call == 1.0:
                prob *= emit_p
            else:
                prob *= 1 - emit_p
        return prob
 
    def _call_row(self, row):
        calls = numpy.zeros(self.M)
        i = 0
        for name in common.DATA_SETS:
            if row[name] >= self.thresholds[name]:
                calls[i] = 1.0
            i += 1
        return calls
 
def find_empirical_means(tab):
    means = {}
    num_windows = len(tab.cols.start)
    for name in common.DATA_SETS:
        means[name] = sum(row[name] for row in tab)/num_windows
    return means

def find_threshold_vals(means):
    thresholds = {}
    for name in common.DATA_SETS:
        for i in range(1,10):
            mean = means[name]
            if stats.poisson.pmf(i, mean) < 1E-4: 
                thresholds[name] = i
                break
    return thresholds

def init_hmm(tab, K):
    means = find_empirical_means(tab)
    thresholds = find_threshold_vals(means)    
    return HMM(tab, means, thresholds, K)
