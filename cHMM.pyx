import ctypes
import numpy as np
cimport numpy as np

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

    def random_init(self):
        self.trans_probs = np.random.random((self.K, self.K))
        self.start_probs = np.random.random(self.K) 
        self.emission_probs = np.random.random((self.K, self.M))

    def prob(self, obs):
        pass

    def _count_transitions(self, obs):
        A = np.zeros((self.K, self.K))
        pass

    def viterbi_decode(self, obs):
        path = np.zeros(len(obs))
        pass

    def all_obs_prob(self, obs):
        probs = np.zeros(self.K)
        for s in range(0, self.K):
            probs[s] = obs_prob(self.emission_probs, s, obs)
        return probs

 

    # AFter "Numerically Stable Hidden markov Model Implementation" by Tobias Mann
cdef forward(hmm, np.ndarray[np.float64_t, ndim=2] obs):
    cdef int K = hmm.K
    cdef int T = obs.shape[0]
    cdef np.ndarray[dtype_t, ndim=2] alpha = np.zeros((T, K), dtype=np.float64)
    cdef np.ndarray[dtype_t, ndim=1] scale = np.zeros(T, dtype=np.float64)
    cdef np.ndarray[dtype_t, ndim=2] trans_probs = hmm.trans_probs
    cdef np.ndarray[dtype_t, ndim=1] start_probs = hmm.start_probs
    cdef Py_ssize_t t, i, j
    alpha[0] = hmm.start_probs * hmm.all_obs_prob(obs[0])
    scale[0] = sum(alpha[0])
    alpha[0] /= scale[0]
#        alpha[0, :] = numpy.ones(self.K)
    for t from 1 <= t < T: 
        obs_probs = hmm.all_obs_prob(obs[t])
        alpha[t] = np.array([sum(alpha[t-1] * trans_probs[:,j]) for j in range(0, K)]) * obs_probs
        s = sum(alpha[t, :])
        scale[t] = s
        alpha[t] = alpha[t, :]/s
    print alpha[T-1] 
    return (alpha, scale)
 
cdef backward(hmm, np.ndarray[dtype_t, ndim=2] obs, np.ndarray[dtype_t, ndim=1] scale):
    cdef int T = obs.shape[0]
    cdef int K = hmm.K
    cdef np.ndarray[np.float64_t, ndim=2] trans_probs = hmm.trans_probs
    cdef np.ndarray[np.float64_t, ndim=2] beta = np.zeros((T, K), dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] obs_probs = hmm.all_obs_prob(obs[T-1])
    cdef Py_ssize_t j, t, i
#        for j from 0 <= j < self.K:
#            for i from 0 <= i < self.K:
    beta[T - 1, i] = np.array([sum(trans_probs[j,:] * obs_probs) for j in range(0, K)]) / scale[T - 1]
        #beta[T - 1, :] = self.start_probs
    for t from T-2 >= t > -1 by -1:
        obs_probs = hmm.all_obs_prob(obs[t+1])
        for j from 0 <= j < K:
            a = 0
            for i from 0 <= i < K:
                a += beta[t+1][i]  
            beta[t, j] = np.array([sum(beta[t+1] * obs_probs * trans_probs[j,:]) for j in range(0, K)]) / scale[t]
        #return sum(self.starts_probs * self.all_obs_prob(obs[0]) * beta[0])
    print beta[0]
    return beta

     # Returns an array of the likelihood of emitting obs for all hidden states  

        # Determines probaility distribution over set of observations 
cdef obs_prob(np.ndarray[dtype_t, ndim=2] emission_probs, int state, np.ndarray[dtype_t, ndim=1] obs):
    cdef np.ndarray[dtype_t, ndim=1] call
    cdef Py_ssize_t i
    cdef dtype_t prob
    prob = 1.0
    for i from 0 <= i < len(obs):
        call = obs[i] # must be either 0 or 1
        emit_p = emission_probs[state, i] # prob of mark i in state
        if call == 1.0:
            prob *= emit_p
        else:
            prob *= 1 - emit_p
    return prob
