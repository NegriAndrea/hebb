#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def weighted_median(values, weights):
    i = np.argsort(values)
    c = np.cumsum(weights[i])
    return values[i[np.searchsorted(c, 0.5 * c[-1])]]

def bootstrap(data, niterations):
    NMassRank = data.shape[0]
    medians = np.zeros((NMassRank, niterations), dtype=data.dtype)
    for j in range(niterations):
        indexes=np.random.randint(0, data.shape[1], size=data.shape[1])
        Mtest=data[:,indexes]
        for k in range(NMassRank):
            medians[k,j] = np.median(Mtest[k,:])

    return medians

def bbootstrap(data, niterations):
    NMassRank = data.shape[0]
    medians = np.zeros((NMassRank, niterations), dtype=data.dtype)
    ones = np.ones(data.shape[1], dtype=np.float32)
    for j in range(niterations):
        weights = np.random.dirichlet(ones)
        for k in range(NMassRank):
            medians[k,j] = weighted_median(data[k,:], weights)

    return medians
