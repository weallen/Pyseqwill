import ctypes
import numpy as np
cimport numpy as np
cdef extern from "math.h":
    cdef extern float exp(float x)
    cdef extern float log(float x)
    cdef extern float sqrt(float x)
    cdef extern float pow(float x, float exponent)

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

    def set_obs(self, obs):
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
        self.emission_probs = np.ones((self.K, self.M))

    def prob(self, obs):
        pass

    def forward(self):
        return chmm_forward(self)
    
    def backward(self, scale):
        return chmm_backward(self, scale)

    def _count_transitions(self, obs):
        A = np.zeros((self.K, self.K))
        pass

    def decode(self):
        return chmm_decode(self)
    

cdef chmm_forward_backward(hmm):
    cdef int T = hmm.all_obs_probs.shape[0]
    cdef int K = hmm.K
    cdef np.ndarray[dtype_t, ndim=1] scale = np.zeros(K, dtype=np.float64)
    cdef np.ndarray[dtype_t, ndim=2] alpha = np.zeros((T, K), dtype=np.float64)
    cdef np.ndarray[dtype_t, ndim=2] beta = np.zeros((T, K), dtype=np.float64)
    cdef np.ndarray[dtype_t, ndim=2] gamma = np.zeros((T, K), dtype=np.float64)
    cdef np.ndarray[dtype_t, ndim=2] eta = np.zeros((T, K), dtype=np.float64)
    int iters = 0
    float prev_loglik = -999999999.0
    float loglik = chmm_forward(hmm, alpha, scale)
    while  iters < common.MAX_ITERS and loglik > prev_loglik:
        loglik = chmm_forward(hmm, alpha, scale)
        chmm_backward(hmm, beta, scale)
        chmm_state_prob(hmm, gamma, alpha, beta) 
        chmm_reestimate_parameters(hmm, eta, alpha, beta)
        prev_loglik = loglik  
        iters += 1

cdef chmm_decode(hmm):
    return 1


cdef chmm_reestimate_parameters(hmm, np.ndarray[dtype_t, ndim=2] eta, 
                np.ndarray[dtype_t, ndim=2] alpha, 
                np.ndarray[dtype_t, ndim=2] beta):    
    pass
  
cdef chmm_state_prob(hmm, np.ndarray[dtype_t, ndim=2] gamma,
                np.ndarray[dtype_t, ndim=2] alpha, 
                np.ndarray[dtype_t, ndim=2] beta):
    cdef int K = hmm.K
    cdef int T = hmm.all_obs_probs.shape[0]
    cdef np.ndarray[dtype_t, ndim=2] gamma = np.zeros((T, K), dtype=np.float64)
    int t, norm, i
    for t from 0 <= t < T:
        for i from 0 <= i < K:
            gamma[t, i] = alpha[t, i] * beta[t, i]
#        chmm_normalize_in_place(gamma, t, K)
   
# AFter "Numerically Stable Hidden markov Model Implementation" by Tobias Mann
cdef chmm_forward(hmm, np.ndarray[dtype_t, ndim=2] alpha, np.ndarray[dtype_t, ndim=1] scale):
    cdef int K = hmm.K
    cdef int T = hmm.all_obs_probs.shape[0]
    cdef np.ndarray[dtype_t, ndim=2] trans_probs = hmm.trans_probs
    cdef np.ndarray[dtype_t, ndim=1] start_probs = hmm.start_probs
    cdef np.ndarray[dtype_t, ndim=2] all_obs_probs = hmm.all_obs_probs
    cdef int t, i, j, s
    cdef dtype_t a
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
    float loglik = 0;
    for t from 0 <= t < T: 
        loglik += log(scale[t])
    return -loglik
 
cdef chmm_backward(hmm, np.ndarray[dtype_t, ndim=2] beta):
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
cdef chmm_all_obs_prob(hmm, np.ndarray[long, ndim=2] obs):
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
cdef chmm_obs_prob(np.ndarray[dtype_t, ndim=2] emission_probs, int state, np.ndarray[long, ndim=1] obs):
    cdef long call
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
 
