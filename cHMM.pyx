import ctypes
import numpy as np
cimport numpy as np

from scipy import stats
import numpy
import common
import table
import math
import common
cdef extern from "math.h":
    cdef float log(float x)


ctypedef np.float64_t dtype_t

# K = number of hidden states
# M = number of input marks
class HMM:
    def __init__(self, M, K):
        self.M = M
        self.K = K
        # prob that will transition from state i to state j
        self.trans_probs = np.zeros((self.K, self.K), dtype=np.float64)
        # prob that in state k has mark m
        self.emission_probs = np.zeros((self.K, self.M), dtype=np.float64) 
        self.start_probs = np.zeros(self.K, dtype=np.float64)
        self.all_obs_probs = None
        self.obs = None

    def set_obs(self, obs):
        print "Finding empirical means"
        means = find_empirical_means(obs)
        print means
        print "Finding thresholds"
        thresholds = find_threshold_vals(means)
        print thresholds
        print "Findng obs seq"
        self.obs = chr_obs_seq(obs, thresholds)
        print "Finding obs prob"
        self.all_obs_probs = chmm_all_obs_prob(self, obs)

    def random_init(self):
        self.trans_probs = np.random.random((self.K, self.K))
        for i in range(self.K):
            norm = sum(self.trans_probs[i])
            for j in range(self.K):
                self.trans_probs[i, j] /= norm
        self.start_probs = np.random.random(self.K) 
        norm = sum(self.start_probs)
        for i in range(self.K):
            self.start_probs[i] /= norm
        self.emission_probs = np.ones((self.K, self.M))/self.M

    def train(self):
        chmm_forward_backward(self)

    
cdef inline int is_over_thresh(val, thresh):
    return (1 if val >= thresh else 0)

cdef call_row(np.ndarray[double, ndim=2] chr_cov, int row_idx, 
                np.ndarray[dtype_t, ndim=1] thresholds, int row_len):
    cdef np.ndarray[double, ndim=1] calls = np.zeros(row_len)
    for i in range(row_len):
        calls[i] = is_over_thresh(chr_cov[row_idx, i], thresholds[i]) 
    return calls

# chr_cov is 2D (num datasets x num windows) matrix of window counts
cdef chr_obs_seq(np.ndarray[double, ndim=2] chr_cov, np.ndarray[dtype_t, ndim=1] thresholds):
    cdef int num_rows = len(common.DATA_SETS)
    cdef int num_obs = chr_cov.shape[0]
    cdef int i, j
    cdef np.ndarray[double, ndim=1] currrow
    cdef obs = numpy.zeros((num_obs, num_rows), dtype=np.int32)
    for i from 0 <= i < num_obs:
        currrow = call_row(chr_cov, i, thresholds, num_rows)  
        for j from 0 <= j < num_rows:
            obs[i,j] = currrow[j]
    return obs

def find_empirical_means(obs):
    means = np.zeros(len(common.DATA_SETS))
    for i in range(len(common.DATA_SETS)):
        means[i] = round(np.mean(obs[:, i]))
    return means

def find_threshold_vals(means):
    thresholds = np.zeros(means.shape[0])
    for d in range(len(common.DATA_SETS)):
        for i in range(1,15):
            mean = means[d]
            if stats.poisson.pmf(i, mean) < 1E-4: 
                thresholds[d] = i
                break
    return thresholds
   
def chmm_forward_backward(hmm):
    T = hmm.all_obs_probs.shape[0]
    K = hmm.K
    scale = np.zeros(T, dtype=np.float64)
    alpha = np.zeros((T, K), dtype=np.float64)
    beta = np.zeros((T, K), dtype=np.float64)
    gamma = np.zeros((T, K), dtype=np.float64)
    eta = np.zeros((T, K, K), dtype=np.float64)
    prev_loglik = -999999999.0
    loglik = chmm_forward(hmm, alpha, scale)
    print loglik
    iters = 0
    while iters < common.MAX_ITERS and loglik > prev_loglik:
        if iters > 0:
            loglik = chmm_forward(hmm, alpha, scale)
        print loglik
        chmm_backward(hmm, beta, scale)
        chmm_state_prob(hmm, eta, gamma, alpha, beta) 
        chmm_reestimate_parameters(hmm, eta, gamma, alpha, beta)
        prev_loglik = loglik  
        iters += 1

cdef chmm_state_prob(hmm, np.ndarray[dtype_t, ndim=3] eta,
                np.ndarray[dtype_t, ndim=2] gamma, 
                np.ndarray[dtype_t, ndim=2] alpha, 
                np.ndarray[dtype_t, ndim=2] beta ):    
    cdef int K = hmm.K
    cdef int T = hmm.all_obs_probs.shape[0]
    cdef np.ndarray[dtype_t, ndim=2] trans_probs = hmm.trans_probs
    cdef np.ndarray[dtype_t, ndim=1] start_probs = hmm.start_probs
    cdef np.ndarray[dtype_t, ndim=2] all_obs_probs = hmm.all_obs_probs
    cdef int t, i, j
    cdef float denom
    for t from 0 <= t < T-1:
        denom = 0.0
        for i from 0 <= i < K:
            for j from 0 <= j < K:
                denom += alpha[t, i] * trans_probs[i, j] * all_obs_probs[t + 1, j] * beta[t + 1, j]
        for i from 0 <= i < K:
            gamma[t, i] = 0.0
            for j from 0 <= j < K:
                eta[t, i, j] = (alpha[t, i] * trans_probs[i, j] * all_obs_probs[t + 1, j] * beta[t + 1, j])/denom
                gamma[t, i] += eta[t, i, j]

cdef chmm_reestimate_parameters(hmm, np.ndarray[dtype_t, ndim=3] eta,
                np.ndarray[dtype_t, ndim=2] gamma,
                np.ndarray[dtype_t, ndim=2] alpha, 
                np.ndarray[dtype_t, ndim=2] beta):
    cdef int K = hmm.K
    cdef int M = hmm.M
    cdef int T = hmm.all_obs_probs.shape[0]
    cdef np.ndarray[dtype_t, ndim=2] trans_probs = hmm.trans_probs
    cdef np.ndarray[dtype_t, ndim=1] start_probs = hmm.start_probs
    cdef np.ndarray[dtype_t, ndim=2] all_obs_probs = hmm.all_obs_probs 
    cdef np.ndarray[int, ndim=2] obs = hmm.obs
    cdef np.ndarray[dtype_t, ndim=2] emission_probs = hmm.emission_probs
    cdef int t, norm, i, j
    cdef float numer, denom
    # reestimate start_probs
    for i from 0 <= i < K:
        start_probs[i] = gamma[0, i]
    # reestimate trans_probs
    for i from 0 <= i < K:
        for j from 0 <= j < K:
            numer = 0.0
            denom = 0.0
            for t from 0 <= t < T-1:
                numer += eta[t, i, j]
                denom += gamma[t, i]
            trans_probs[i, j] = numer / denom
    # reestimate emission_probs
    for i from 0 <= i < K:
        for j from 0 <= j < M:
            numer = 0.0
            denom = 0.0
            for t from 0 <= t < T - 1:
                if (obs[t, j] == 1):
                    numer += gamma[t, i]
                denom += gamma[t, i]
            emission_probs[i, j] = numer / denom    

cdef chmm_forward(hmm, np.ndarray[dtype_t, ndim=2] alpha, np.ndarray[dtype_t, ndim=1] scale):
    cdef int K = hmm.K
    cdef int T = hmm.all_obs_probs.shape[0]
    cdef np.ndarray[dtype_t, ndim=2] trans_probs = hmm.trans_probs
    cdef np.ndarray[dtype_t, ndim=1] start_probs = hmm.start_probs
    cdef np.ndarray[dtype_t, ndim=2] all_obs_probs = hmm.all_obs_probs
    cdef int t, i, j, s
    cdef float loglik = 0.0
    for i from 0 <= i < K:
        alpha[0, i] = start_probs[i] * all_obs_probs[0, i] 
    scale[0] = chmm_normalize_in_place(alpha, 0, K)
    for t from 1 <= t < T: 
        for j from 0 <= j < K:
            a = 0
            for i from 0 <= i < K:
                a += alpha[t-1, i] * trans_probs[i, j]
            alpha[t, j] = a * all_obs_probs[t, j]
#        alpha[t] = np.array([sum(alpha[t-1] * trans_probs[:,j]) for j in range(0, K)]) * obs_probs
        scale[t] = chmm_normalize_in_place(alpha, t, K)
    for t from 0 <= t < T: 
        loglik += log(scale[t])
    return -loglik
 
cdef chmm_backward(hmm, np.ndarray[dtype_t, ndim=2] beta, np.ndarray[dtype_t, ndim=1] scale):
    cdef int T = hmm.all_obs_probs.shape[0]
    cdef int K = hmm.K
    cdef np.ndarray[np.float64_t, ndim=2] trans_probs = hmm.trans_probs
    cdef np.ndarray[np.float64_t, ndim=2] all_obs_probs = hmm.all_obs_probs
    cdef int j, t, i
#    beta[T - 1, i] = np.array([sum(trans_probs[j,:] * obs_probs) for j in range(0, K)]) / scale[T - 1]
    for i from 0 <= i < K: 
        beta[T - 1, i] = 1/scale[T-1]
    for t in range(T-2, -1, -1):
        for i from 0 <= i < K:
            b = 0
            for j from 0 <= j < K:
                b += trans_probs[i,j] * all_obs_probs[t+1, j] * beta[t+1, j]
#            beta[t] = np.array([sum(beta[t+1] * obs_probs * trans_probs[j,:]) for j in range(0, K)]) / scale[t]
            beta[t, i] = b / scale[t]

# Returns an array of the likelihood of emitting each obs for all hidden states  
cdef chmm_all_obs_prob(hmm, np.ndarray[double, ndim=2] obs):
    cdef int T = obs.shape[0]
    cdef int K = hmm.K
    cdef int i, j
    cdef np.ndarray[dtype_t, ndim=2] probs = np.zeros((T, K), dtype=np.float64)
    cdef np.ndarray[dtype_t, ndim=2] emission_probs = hmm.emission_probs
    for i from 0 <= i < T:
        for j from 0 <= j < K:
            probs[i, j] = chmm_obs_prob(emission_probs, j, obs[i])
    return probs


# Determines probaility distribution over set of observations 
cdef chmm_obs_prob(np.ndarray[dtype_t, ndim=2] emission_probs, int state, np.ndarray[double, ndim=1] obs):
    cdef double call
    cdef int i
    cdef float prob = 1.0
    for i from 0 <= i < len(obs):
        call = obs[i] # must be either 0 or 1
        emit_p = emission_probs[state, i] # prob of mark i in state
        if call == 1:
            prob *= emit_p
        else:
            prob *= 1 - emit_p
    return prob

cdef chmm_normalize_in_place(np.ndarray[dtype_t, ndim=2] A, int i, int N):
    cdef int n
    cdef float s = 0
    for n from 0 <= n < N:
        s += A[i, n]
    for n from 0 <= n < N:
        A[i, n] /= s
    return s
 
