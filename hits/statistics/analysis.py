import re, csv, numpy, math, scipy
import scipy.optimize
import scipy.stats
import numpy as np
from itertools import repeat
from .misc_para import parmap

#This is essentially used as a closure
class LogLikelihood(object):    
    def __init__(self, theta, beta, x, y):
        self.theta = theta
        self.beta = beta
        self.x = x
        self.y = y

    def __call__(self, m):
        r1 = 1/self.theta
        r2 = self.beta * r1
        p = 1/(1+(self.theta*m))
        l1 = scipy.stats.nbinom.logpmf(self.x, r1, p)
        l2 = scipy.stats.nbinom.logpmf(self.y, r2, p)
        return -(l1 + l2)

def _compute_single(inputs):
    '''
    Uses given theta/beta parameters to find the mu parameter most likely to produce the observed frequencies under a NB distribution and returns the likelihood
    '''
    theta, beta, x, y = inputs
    ll = LogLikelihood(theta, beta, x, y)
    bnds = (1e-8, 1e6)
    return scipy.optimize.minimize_scalar(ll, bounds = bnds, method = 'bounded').fun #minimized objective (minimized by mu)

def _compute_single_m(inputs):
    '''
    Identify the most likely mean parameter to generate the observed counts for a specific variant in the 2 libraries given the theta/beta parameters
    '''
    theta, beta, x, y = inputs
    ll = LogLikelihood(theta, beta, x, y)
    bnds = (1e-8, 1e6)
    m = scipy.optimize.minimize_scalar(ll, bounds = bnds, method = 'bounded').x #mu parameter giving minimized objective
    return m

def compute_m(theta, beta, sample1, sample2):
    '''
    Identify the most likely mean parameters to generate the observed counts for all variants in the 2 libraries given the theta/beta parameters
    '''
    return parmap(_compute_single_m, zip(repeat(theta), repeat(beta), sample1, sample2))

## Function for computing p-value using negative binomial distribution ##
def _compute_single_pv(inputs):
    '''
    Find P-value with one-tailed exact test based on fitted NB distribution conditioned on total observed count in both libraries 
    '''
    # One-tailed test
    theta, beta, x, y, m, greater = inputs
    r1 = 1/theta
    r2 = beta * r1
    p = 1/(1+(theta*m))
    l1 = scipy.stats.nbinom.logpmf(x, r1, p)
    l2 = scipy.stats.nbinom.logpmf(y, r2, p)
    ll = l1 + l2
    sum_ll = ll 
    sum_small = ll
    s = x + y 
    for j in range(0, int(s)):
        if (j == y):
            continue
        l1 = scipy.stats.nbinom.logpmf(s-j, r1, p)
        l2 = scipy.stats.nbinom.logpmf(j, r2, p)
        l = l1 + l2
        sum_ll = numpy.logaddexp(sum_ll, l)
        #if l < ll:
        #greater is boolean flag for which single tailed test
        if greater and j>y:
            sum_small = numpy.logaddexp(sum_small, l)
        elif not greater and j<y:
            sum_small = numpy.logaddexp(sum_small, l)
    log10_pv = (sum_small - sum_ll)/numpy.log(10)

    return log10_pv 

def compute_pv(theta, beta, sample1, sample2, m, greater=True):
    return parmap(_compute_single_pv, zip(repeat(theta), repeat(beta), sample1, sample2, m, repeat(greater)))

def BH_correction_log(pv, alpha):
    '''
    Adjust set of P-values for multiple hypothesis testing
    '''
    #Benjamani Hochberg correction for log10 pvalues
    #sort pvalues and store originial index
    pv_sort = numpy.sort(pv)
    pv_index = numpy.argsort(pv)

    #Calculate BH critical values
    BH_critical = [(float(i)/len(pv_sort))*alpha for i in range(1, len(pv_sort)+1)]

    #Count significant values 
    largest_index = []
    for i in range(0, len(pv_sort)):
        if pv_sort[i] <= numpy.log10(BH_critical[i]):
            largest_index.append(i)

    if len(largest_index) >= 1:
        nsig = max(largest_index)+1
    else:
        nsig = 0

    #Adjust pvalues and reindex to original index
    pv_adj = []
    for i in range(1, len(pv_sort)+1):
        if i < len(pv_sort):
            if (pv_sort[i-1]+numpy.log10(float(len(pv_sort))/i)) <= (pv_sort[i]+numpy.log10(float(len(pv_sort))/(i+1))):
                pv_adj.append(pv_sort[i-1]+numpy.log10(float(len(pv_sort))/i))
            if (pv_sort[i-1]+numpy.log10(float(len(pv_sort))/i)) > (pv_sort[i]+numpy.log10(float(len(pv_sort))/(i+1))):
                pv_adj.append((pv_sort[i]+numpy.log10(float(len(pv_sort))/(i+1))))
        if i == len(pv_sort):
            pv_adj.append(pv_sort[i-1]+numpy.log10(float(len(pv_sort))/i))

    pv_adj_reindex = numpy.zeros(len(pv_index))
    for i in range(0, len(pv_index)):
        if float(pv[pv_index[i]]) == 0.0:
            pv_adj_reindex[pv_index[i]] = 0.0
        if float(pv[pv_index[i]]) == 1.0:
            pv_adj_reindex[pv_index[i]] = 1.0
        if float(pv[pv_index[i]]) != 0.0 and float(pv[pv_index[i]]) != 1.0: 
            pv_adj_reindex[pv_index[i]] = pv_adj[i]

    return pv_adj_reindex

class FitParams(object):
    def __init__(self, data):
        self.sample1 = data[0]
        self.sample2 = data[1]

    def _compute0(self, theta, beta):
        like = 0
        for x,y in zip(self.sample1, self.sample2):
            ll = LogLikelihood(theta, beta, x, y)

            #initial guess, matching R optimize function 
            bnds = (min((x,y)), max((x,y)))
            opt = scipy.optimize.minimize_scalar(ll, bounds = bnds, method = 'bounded') 
            like = like + opt.fun

        return like

    def _compute(self, theta, beta):
        '''
        Sums likelihood over the variants.
        '''
        return np.sum(parmap(_compute_single, \
                zip(repeat(theta), repeat(beta), self.sample1, self.sample2)))

    def __call__(self, params):
        theta, beta = params
        return self._compute(theta, beta)

class FitParamsBeta(FitParams):
    def __init__(self, data, theta):
        self.theta = theta
        super(FitParamsBeta, self).__init__(data)

    def __call__(self, params):
        beta = params
        return self._compute(self.theta, beta)

class FitParamsTheta(FitParams):
    def __init__(self, data, beta):
        self.beta = beta
        super(FitParamsTheta, self).__init__(data)

    def __call__(self, params):
        theta = params
        return self._compute(theta, self.beta)

def fit_background(data, theta=None, beta=None):
    '''
    Finds the theta/beta parameters most consistent with the observed data
    '''
    theta_bounds = (1e-8, 10.)
    beta_bounds = (0.001, 15.)

    if theta is None and beta is None:
        fit_fn = FitParams(data)
        bnds = (theta_bounds, beta_bounds)
        guess = (1e-8, 1.0)
        params = scipy.optimize.minimize(fit_fn, guess, bounds=bnds).x 
        theta, beta = params
    elif theta is None:
        fit_fn = FitParamsTheta(data, beta)
        guess = 1e-8
        theta = scipy.optimize.minimize(fit_fn, guess, bounds=(theta_bounds,)).x[0]
    elif beta is None:
        fit_fn = FitParamsBeta(data, theta)
        guess = 1.
        beta = scipy.optimize.minimize(fit_fn, guess, bounds=(beta_bounds,)).x[0] 

    ### Output the estimated values for beta and gamma
    return theta, beta

def analyze_samples(sample1, sample2, theta, beta):
    '''
    Finds the statistical significance of the provided samples given the global theta/beta parameters
    '''
    m = compute_m(theta, beta, sample1, sample2)
    pv_greater = compute_pv(theta, beta, sample1, sample2, m)
    pv2_greater = BH_correction_log(pv_greater, alpha=0.001)

    print("There are %d significant hits" % np.sum(pv2_greater<-3))

    #observation_threshold = 0
    for observation_threshold in range(100):
        m = compute_m(theta, beta, [0]*len(sample1), [observation_threshold]*len(sample1))
        pv = compute_pv(theta, beta, [0]*len(sample1), [observation_threshold]*len(sample1), m)
        pv2 = BH_correction_log(pv, alpha=0.001)
        #print(pv2[0])
        if pv2[0] < -3: break

    #pv_lesser = compute_pv(theta, beta, sample1, sample2, m, greater=False)
    #pv2_lesser = BH_correction_log(pv_lesser, alpha=0.001)

    #pv2 = np.minimum(pv2_greater, pv2_lesser)
    pv2 = pv2_greater

    index = np.where(pv2_greater<-3)[0]
    log_fc = np.log2(((sample2+1)/(sample1+1).astype(float)))
    #observed_mask = (sample2>0) | (sample1>0)
    observed_mask = (sample2 + sample1) > observation_threshold
    print("Coverage = %.3f of variants with at least %d total reads which is minimum required for statistical significance"%(np.sum(observed_mask)/float(len(observed_mask)), observation_threshold))

    unobserved_mask = np.invert(observed_mask)
    observed = np.where(observed_mask)[0]+1
    unobserved = np.where(unobserved_mask)[0]+1
    pv2[unobserved-1] = 10
    return log_fc, observed_mask, pv2, observed, unobserved
