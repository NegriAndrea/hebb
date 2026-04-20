#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import numpy.typing as npt
from astropy.table import Table
from time import perf_counter
from .logger_mine import loggerH, loggerN

def weighted_median(values: npt.NDArray, weights: npt.NDArray) -> npt.NDArray:
    i = np.argsort(values)
    c = np.cumsum(weights[i])
    return values[i[np.searchsorted(c, 0.5 * c[-1])]]

def bootstrap(data: npt.NDArray, niterations: int) -> npt.NDArray:
    NMassRank = data.shape[0]
    medians = np.zeros((NMassRank, niterations), dtype=data.dtype)
    for j in range(niterations):
        indexes=np.random.randint(0, data.shape[1], size=data.shape[1])
        Mtest=data[:,indexes]
        for k in range(NMassRank):
            medians[k,j] = np.median(Mtest[k,:])

    return medians

# keep this isolated from bootstrap_driver in case I want to use numba
# for now quantiles with weights is not numba supported
def bbootstrap(data: npt.NDArray, quantiles: tuple, niterations: int) -> npt.NDArray:
    """
    Perform a bayesian bootstrap on the quantiles of the data.

    """

    NMassRank = data.shape[0]
    Q = np.zeros((NMassRank, niterations, 3), dtype=data.dtype)
    ones = np.ones(data.shape[1], dtype=np.float32)

    for j in range(niterations):
        weights = np.random.dirichlet(ones)
        for k in range(NMassRank):
            Q[k,j,:] = np.quantile(data[k,:], quantiles, weights = weights,
                                   method="inverted_cdf")

    return Q

def bootstrap_driver(Mmax: npt.NDArray,
                     niterations: int,
                     quantiles: tuple =  (0.16, 0.5, 0.86)) -> Table:
    """
    Perform a bayesian bootstrap of the Mmax quantiles, separately for each
    mass rank.

    Parameters
    ----------
    Mmax : np.ndarray of shape (NMassRanks, Nboxes)
        A 2D array of floats representing the maximum halo masses
        found in each bootstrap iteration. Units should be log10(M_sun).

    niterations: int
        Number of iterations of the bayesian bootstrap

    quantiles: tuple, optional
        The quantiles used to compute the bootstrap, default (0.16, 0.5, 0.86)


    Returns
    -------
    Astropy table of variance of the log quantiles, for every mass rank.

    """
    loggerH(f"BAYESIAN BOOTSTRAP: beginning on a sample of {Mmax.shape[1]}"
            " boxes")
    t0 = perf_counter()
    Mmax_float32 = (10.**Mmax.astype(np.float64)).astype(np.float32)
    quantilesValues = bbootstrap(Mmax_float32, quantiles, niterations)

    out = np.sqrt(np.var(np.log10(quantilesValues), axis = 1))

    loggerN(f"BAYESIAN BOOTSTRAP: took {perf_counter()-t0:.1f} s")
    return out
