"""Implement uncertainty metrics from 
J. T. Horwood, J. M. Aristoff, N. Singh, A. B. Poore, and M. D. Hejduk,
"Beyond covariance realism: a new metric for uncertainty realism"
in Proceedings of the SPIE, Signal and Data Processing of Small Targets, 
Vol. 9092, Baltimore, MD, May 2014

Basic idea is to compute chi^2 given residuals and a covariance matrix.
Then test the distribution of chi^2 against a chi^2 distribution, whether
by binning up the frequencies of occurrence and comparing against predictions
(Pearson's Chi-Squared Goodness-of-Fit test) or by whole-distribution tests
(KS test, Cramer-von Mises, Anderson-Darling).

I'll start with a 5 bin Pearson's Chi-Squared test and a Cramer-von Mises test.
"""

import numpy as np
import scipy.stats


def chi2(dx, covar):
    """Compute chi-squared statistic given residuals and covariance matrices.
    
    Calculates chi^2 = dx^T * C^-1 * dx for each residual vector and
    corresponding covariance matrix.
    
    Parameters
    ----------
    dx : array_like(n, m) or (m,), float
        n residual vectors of dimension m, or single residual vector
    covar : array_like(n, m, m) or (m, m), float
        n covariance matrices (m by m), or single covariance matrix
        
    Returns
    -------
    chi2 : array_like(n) or float
        Array of chi-squared values, or single value if scalar inputs
    """
    if len(covar.shape) == 2 and len(dx.shape) == 1:
        scalar = True
        covar = covar[None, ...]
        dx = dx[None, ...]
    elif len(covar.shape) == 3 and len(dx.shape) == 2:
        scalar = False
    else:
        raise ValueError('unrecognized shapes of dx, covar')
    cinv = np.array([np.linalg.inv(c) for c in covar])
    chi2 = np.einsum('ij,ijk,ik->i', dx, cinv, dx)
    if scalar:
        chi2 = chi2[0]
    return chi2


def cvm_chi2_test(chi2, ndof, alpha=False):
    """Perform Cramér-von Mises test for chi-squared distribution fit.
    
    Tests whether the provided chi-squared values follow a chi-squared
    distribution with the specified degrees of freedom using the
    Cramér-von Mises criterion, a goodness-of-fit test based on the
    empirical distribution function.
    
    Parameters
    ----------
    chi2 : array_like(n), float
        Chi-squared values to test
    ndof : int
        Number of degrees of freedom for the reference chi-squared distribution
    alpha : bool, optional
        If True, convert test statistic to approximate p-value
        (probability of obtaining this statistic if the null hypothesis is true).
        Default is False.

    Returns
    -------
    ts : float
        Cramér-von Mises test statistic. Larger values indicate worse fit.
        If alpha=True, returns approximate p-value instead.
        
    Notes
    -----
    The Cramér-von Mises test is more sensitive to differences in the middle
    of the distribution compared to the Kolmogorov-Smirnov test.
    """
    chi2 = np.asarray(chi2)
    chi2 = chi2[np.argsort(chi2)]
    cdf = scipy.stats.chi2.cdf(chi2, ndof)
    n = len(chi2)
    count = np.arange(n)+1
    ts = 1/(12*n) + np.sum(((2*count-1)/(2*n) - cdf)**2)
    if alpha:
        return cvm_to_alpha(ts, n)
    else:
        return ts


def cvm_to_alpha(ts, n):
    """Convert Cramér-von Mises test statistic to approximate p-value.
    
    Uses the asymptotic critical values from Stephens (1986) with a
    correction for finite sample sizes to estimate the p-value.

    Parameters
    ----------
    ts : float
        Cramér-von Mises test statistic
    n : int
        Number of data points used in computing the test statistic.
        Must be >= 5.

    Returns
    -------
    alpha : float
        Approximate p-value: the probability of obtaining a test statistic
        this large or larger if the data were drawn from the reference
        distribution. Values closer to 1 indicate good fit, values closer
        to 0 indicate poor fit.
        
    Notes
    -----
    Based on Table 4.2 from:
    M. A. Stephens (1986). "Tests Based on EDF Statistics". 
    In D'Agostino, R.B.; Stephens, M.A. (eds.). Goodness-of-Fit Techniques
    
    Values for alpha > 0.25 are mostly intended to be illustrative and
    should be interpreted with caution.
    
    Raises
    ------
    ValueError
        If n < 5
    """
    if n < 5:
        raise ValueError('n must be >= 5')
    ts = (ts - 0.4/n + 0.6/n**2) / (1+1/n)
    # table from M. A. Stephens (1986)
    alpha = np.array([0.25, 0.15, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001])
    critval = np.array([0.209, 0.284, 0.347, 0.461, 0.581, 0.743, 0.869, 
                        1.167])
    # extra values for lower alpha
    alpha = np.concatenate([[0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3], alpha])
    critval = np.concatenate([
            [0.025, 0.046, 0.062, 0.078, 0.097, 0.118, 0.146, 0.184], critval])
    out = np.interp(ts, critval, alpha)
    return out


def pearsons_chi(chi2, ndof, nbin):
    """Perform binned Pearson's chi-squared goodness-of-fit test.
    
    Bins the chi-squared values and compares observed bin counts to expected
    uniform distribution across quantiles of the reference chi-squared
    distribution. Returns normalized residuals for each bin.
    
    Parameters
    ----------
    chi2 : array_like(n), float
        Chi-squared values to test
    ndof : int
        Number of degrees of freedom for the reference chi-squared distribution
    nbin : int
        Number of bins to use for the test
        
    Returns
    -------
    normalized_residuals : array_like(nbin), float
        Normalized residuals (observed - expected) / sqrt(expected) for each bin.
        Under the null hypothesis, these should follow approximately a
        standard normal distribution. Large positive values indicate excess
        observations in that bin, large negative values indicate deficit.
        
    Notes
    -----
    The bins are chosen to have equal probability under the reference
    chi-squared distribution (equiprobable binning). This is appropriate
    for testing whether the observed distribution matches the expected
    chi-squared distribution.
    
    The test assumes Poisson statistics in each bin, so the expected count
    in each bin is len(chi2)/nbin with standard deviation sqrt(len(chi2)/nbin).
    """
    edges = np.linspace(0, 1, nbin+1)
    cuts = np.concatenate([[-np.inf], 
                           scipy.stats.chi2.ppf(edges[1:-1], ndof), [np.inf]])
    bin = np.searchsorted(cuts, chi2)-1
    # expected number of counts in each bin is just len(chi2) / nbin
    counts = np.bincount(bin, minlength=nbin)
    expected = len(chi2)/nbin
    # I guess the right thing is a Poisson distribution in each bin.  But I don't
    # really care about the total probability for this particular test;
    # really, we just want to be able to say that there were too few far and
    # too many close observations, IMO.  So let's provide a bunch of chi.
    return (counts-expected)/np.sqrt(expected)