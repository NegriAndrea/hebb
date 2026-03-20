#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def pretty_print(Mmax: np.ndarray) -> None:
    """
    Pretty print a summary of the bootstrap results and quantile table.

    This function takes the extremal mass distribution and computes standard
    quantiles to provide a human-readable summary of the hebb analysis.

    Parameters
    ----------
    Mmax : np.ndarray of shape (NMassRanks, Nboxes)
        A 2D array of floats representing the maximum halo masses
        found in each bootstrap iteration. Units should be log10(M_sun).

    Returns
    -------
    None

    """
    from astropy.table import Table
    quantiles = np.array([0.16, 0.5, 0.86,
                          0.05, 0.1, 0.25, .75, 0.9, 0.95])
    quants = np.zeros((Mmax.shape[0], quantiles.size))

    print('')
    print('RESULTS')
    for i in range(Mmax.shape[0]):
        quants[i,:] = np.quantile(10.**Mmax[i,:].astype(np.float64),
                                            quantiles)
        M16, M50, M86 = quants[i,:3]
        print(f"{i} log10(M200) = {np.log10(M50):.2f} _-{(M50-M16)/M50/np.log(10):.2f} ^+{(M86-M50)/M50/np.log(10):.2f}")

    # build a quantile table and use astropy Table to do a pretty print
    quants = np.log10(quants[:,[3,4,5,1,6,7,8]])
    names=[str(q) for q in quantiles[[3,4,5,1,6,7,8]]]

    t2=Table(quants, names=names)

    # add MassRank column
    t2.add_column(np.arange(quants.shape[0], dtype=np.uint8), name='MassRank',
                  index=0)

    # print only 2 decimals for floating point numbers
    for cname in t2.colnames:
        if t2[cname].info.dtype in ['<f4', '<f8']:
            t2[cname].info.format = '6.2f'

    print('')
    print('QUANTILES TABLE')
    print(t2)
