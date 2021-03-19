#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 14:29:21 2021

@author: Jordan
"""

# Some statistics helper functions
from copy import deepcopy
from math import exp, sqrt, pi
def exp_pdf(x, lam):
    if x < 0.0:
        return 0.0
    else:
        return lam * exp(-lam * x)   
def norm_pdf(x, mu, sigma):
    z = (x - mu) / sigma
    return 1.0 / (sigma * sqrt(2.0 * pi)) * exp(-z * z / 2.0)


##
# A class that represents a model of coverage frequencies in a phased diploid
# assembly, which is fit using data from a coverage frequency histogram in the 
# CoverageDistribution.fit() method
##
class CoverageDistribution:
    
    ############################
    # The downstream interface #
    ############################
    
    # Return the probability that an element in the assembly with coverage c is
    # actually absent in the underlying genome
    def probability_erroneous(self, c):
        return self.component_probabilities(c)[0]
    
    # Return the probability that an element in the assembly with coverage c is
    # present in the underlying genome and haploid (indicating a)
    def probability_haploid(self, c):
        return self.component_probabilities(c)[1]
    
    # Return the probability that an element in the assembly with coverage c is
    # actually consists of multiple elements in the underlying genome collapsed
    # (note that it may be that another copy of this element is where the 
    # collapsing actually ocurred)
    def probability_collapsed(self, c):
        return sum(self.component_probabilities(c)[2:])
    
    # Fit and return a CoverageDistribution using a coverage frequency histogram,
    # which should be provided as either a collections.Counter or dict whose keys
    # are coverages and whose values are the number of bases that have that
    # coverage
    @staticmethod
    def fit(spectrum, tol = 1e-6, max_iters = 10000):
        
        # copy the spectrum so we can modify it without breaking expectations
        spectrum = deepcopy(spectrum)
        
        # This is only necessary with KMC3 output
        # zero out max_val to account for binning at the max value
        #spectrum[max(spectrum)] = 0
        
        # get a coarse estimate that we can use as a starting point for EM
        cov_dist = CoverageDistribution.heuristic_fit(spectrum)
        
        for iteration in range(max_iters):
            print("iter {}: cov={}, var={}, scale={}, weights={}".format(iteration, cov_dist.cov_per_ploidy, sqrt(cov_dist.var_per_ploidy),cov_dist.err_scale, cov_dist.mixture_weights[:3]))
            assignment = {}
            for cov in spectrum:
                assignment[cov] = cov_dist.component_probabilities(cov)
        
            new_mixture_weights = []
            
            # error scale is the weighted mean
            numer = 0.0
            denom = 0.0
            weight_numer = 0.0
            weight_denom = 0.0
            for cov in spectrum:
                freq = spectrum[cov]
                prob = assignment[cov][0]
                numer += freq * prob * cov
                denom += freq * prob
                weight_numer += prob * freq
                weight_denom += freq
                
            new_error_scale = numer / denom
            
            # component weight is weighted proportion assigned to the error component
            new_mixture_weights.append(weight_numer / weight_denom)
                    
            # the haploid coverage is the weighted mean
            numer = 0.0
            denom = 0.0
            for ploidy in range(1, cov_dist.max_ploidy() + 1):
                weight_numer = 0.0
                for cov in spectrum:
                    freq = spectrum[cov]
                    prob = assignment[cov][ploidy]
                    numer += freq * prob * cov
                    denom += prob * freq * ploidy
                    weight_numer += prob * freq
                new_mixture_weights.append(weight_numer / weight_denom)
                
            new_cov_per_ploidy = numer / denom
            
            #print("doing std dev")
            # the haploid std deviation is the sqrt of the weighted average
            # square deviation, adjusted for the fact that the variance of
            # ploidy P should be P times the variance of the haploid
            numer = 0.0
            denom = 0.0
            for cov in spectrum:
                freq = spectrum[cov]
                for ploidy in range(1, cov_dist.max_ploidy() + 1):
                    prob = assignment[cov][ploidy]
                    dev = cov - ploidy * new_cov_per_ploidy
                    numer += freq * prob * dev * dev / ploidy
                    denom += freq * prob
                    
            new_var_per_ploidy = numer / denom
            
            # Limiting the error scale to an expected interval only when 
            # we have a reasonable numebr of zero coverages
            if spectrum[0] > spectrum[1]:
                expected_error_scale =  sqrt(new_var_per_ploidy * 2.0 * pi) * (spectrum[round(new_cov_per_ploidy)] / spectrum[0]) * (new_mixture_weights[0] / new_mixture_weights[1])
                epsilon = abs(new_error_scale - expected_error_scale) / expected_error_scale
                if epsilon > 0.1:
                    new_error_scale = expected_error_scale

            # the next iteration's distrubition
            new_cov_dist = CoverageDistribution(new_cov_per_ploidy, new_var_per_ploidy,
                                                new_error_scale, new_mixture_weights)
            
            #print("doing convergence")
            converged = cov_dist.converged(new_cov_dist, tol)
            cov_dist = new_cov_dist
            if converged:
                break
            
        return cov_dist
    
    ##############################
    # Internal/debugging methods #
    ##############################
    
    # TODO: generalize this for truncated distributions at the max value?
    
    @staticmethod
    def heuristic_fit(spectrum):
        
        max_cov = max(spectrum)
        
        cov_at_max_frequency = None
        max_frequency = -1
        has_increased = False
        for cov_value in sorted(spectrum)[1:]:
            if spectrum[cov_value] > spectrum[cov_value - 1]:
                has_increased = True
            if spectrum[cov_value] > max_frequency and has_increased:
                max_frequency = spectrum[cov_value]
                cov_at_max_frequency = cov_value
        
        cov_per_ploidy = cov_at_max_frequency
        var_per_ploidy = cov_per_ploidy
        
        max_ploidy = int(max_cov / cov_per_ploidy)
        while max_ploidy * cov_per_ploidy - 3.0 * sqrt(max_ploidy * cov_per_ploidy) < max_cov:
            max_ploidy += 1
            
        # start with an arbitrary, small scale
        err_scale = 1.0
        
        # the denominator along with pseudocounts
        total_count = sum(spectrum[cov] for cov in spectrum) + max_ploidy + 1
        
        # approximate the mixture weights simplistically by hard-assigning all data
        # points to their nearest component
        mixture_weights = []
        for i in range(max_ploidy + 1):
            count = 1
            for cov in range(int((i - 0.5) * cov_per_ploidy), int((i + 0.5) * cov_per_ploidy)):
                count += spectrum[cov]
            mixture_weights.append(count / total_count)
        
        return CoverageDistribution(cov_per_ploidy, var_per_ploidy, err_scale, mixture_weights)
    
    def __init__(self, cov_per_ploidy, var_per_ploidy, err_scale, mixture_weights):
        self.cov_per_ploidy = cov_per_ploidy
        self.var_per_ploidy = var_per_ploidy
        self.err_scale = err_scale
        self.mixture_weights = mixture_weights
    
    def pdf(self, x):
        return sum(self.likelihood(ploidy, x) for ploidy in range(self.max_ploidy() + 1))
    
    def max_ploidy(self):
        return len(self.mixture_weights) - 1
    
    def likelihood(self, ploidy, cov):
        l = 0.0
        if ploidy == 0:
            # error distribution, switch to rate parameterization
            l = exp_pdf(cov, 1.0 / self.err_scale)
        else:
            l = norm_pdf(cov, self.cov_per_ploidy * ploidy, sqrt(self.var_per_ploidy * ploidy))
        return self.mixture_weights[ploidy] * l
    
    def component_probabilities(self, cov):
        likelihoods = [self.likelihood(ploidy, cov) for ploidy in range(self.max_ploidy() + 1)]
        Z = sum(likelihoods)
        for i in range(len(likelihoods)):
            likelihoods[i] /= Z
        return likelihoods
    
    def converged(self, other, tol):
        if other.cov_per_ploidy != 0.0 and abs(self.cov_per_ploidy / other.cov_per_ploidy - 1.0) > tol:
            return False
        if other.var_per_ploidy != 0.0 and abs(sqrt(self.var_per_ploidy / other.var_per_ploidy) - 1.0) > tol:
            return False
        if other.err_scale != 0.0 and abs(self.err_scale / other.err_scale - 1.0) > tol:
            return False
        max_ploidy_of_interest = 4
        for wt, other_wt, ploidy in zip(self.mixture_weights, other.mixture_weights, range(len(self.mixture_weights))):
            if ploidy > max_ploidy_of_interest:
                break
            if other_wt and abs(wt / other_wt - 1.0) > tol:
                return False
        return True
    
    
if __name__ == "__main__":
    
    from random import expovariate, normalvariate
    from collections import defaultdict
    import matplotlib.pyplot as plt
    
    ###################################
    # Generate some semi-realistic data
    ###################################
    
    # the simulation parameters
    N = 100000
    var_inflation = 1.5
    max_cov = 255
    hap_cov = 25
    err_frac = .01
    hap_frac = .95
    dip_frac = .01
    high_copy_total_frac = .03
    err_scale = 1.0
    high_copy_fracs = [exp(-cov / hap_cov) for cov in range(3 * hap_cov, 2 * max_cov, hap_cov)]
    Z = sum(high_copy_fracs)
    for i in range(len(high_copy_fracs)):
        high_copy_fracs[i] *= high_copy_total_frac / Z
    mixture_fracs = [hap_frac, dip_frac] + high_copy_fracs
    # do the simulation
    counts = defaultdict(int)
    for i in range(round(err_frac * N)):
        x = int(expovariate(err_scale) + 1)
        if x <= max_cov:
            counts[x] += 1
    for j in range(len(mixture_fracs)):
        for i in range(round(mixture_fracs[j] * N)):
            x = max(1, int(round(normalvariate((j + 1) * hap_cov, sqrt((j + 1) * var_inflation * hap_cov)))))
            if x <= max_cov:
                counts[x] += 1
                
    
    ###################################
    # Fit the model
    ###################################
    
    cov_model = CoverageDistribution.fit(counts, tol = 1e-4, max_iters=1000)

    # look at how good the fit is:
    plt.plot(sorted(counts), [counts[f] for f in sorted(counts)])
    X = list(range(1, max_cov))
    Y = [cov_model.pdf(x) for x in X]
    scale_factor = counts[hap_cov] / Y[hap_cov]
    Y = [scale_factor * y for y in Y]
    plt.plot(X, Y)
    plt.ylim(top = 10000)
    
    ###################################
    # Query some probabilities
    ###################################
    
    for cov in (10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
        prob_err = cov_model.probability_erroneous(cov)
        prob_correct = cov_model.probability_haploid(cov)
        prob_collapsed = cov_model.probability_collapsed(cov)
        
        print("k-mer frequency {}: prob error {:.4f}, hap {:.4f}, coll {:.4f},".format(cov, prob_err, prob_correct, prob_collapsed))
