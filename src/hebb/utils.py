#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import numpy.typing as npt
from .logger_mine import loggerH, loggerN

def pretty_print(Mmax: npt.NDArray) -> None:
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



def save_table(Mmax: npt.NDArray,
               fileNrMax: npt.NDArray,
               subNrMax: npt.NDArray,
               args: argparse.Namespace) -> None:
    """
    Save bootstrap samples in a file.

    Parameters
    ----------
    Mmax : np.ndarray of shape (NMassRanks, Nboxes)
        A 2D array of floats representing the maximum halo masses
        found in each bootstrap iteration. Units should be log10(M_sun).

    fileNrMax: np.ndarray of shape (Nboxes,)
        A 1D array of int containing the IDs of the merger tree files from
        where the subhalos can be found

    subNrMax: np.ndarray of shape (Nboxes,)
        A 1D array of int containing the IDs of the subhalos  in merger tree
        files from where the subhalos can be found

    args: argparse of command line interface
       The arguments of the command line interface, to be saved in the table
       header

    Returns
    -------
    None

    """
    import astropy.units as u
    from astropy.table import Table

    # compute rank in mass and boxID
    massRank = np.repeat(np.arange(Mmax.shape[0]), Mmax.shape[1])
    boxID = np.tile(np.arange(Mmax.shape[1]), Mmax.shape[0])

    # alter shape instead to use np.reshape to force an error
    # in case the code attempts to copy the arrays (should not happen)
    fileNrMax.shape = fileNrMax.size
    subNrMax.shape = subNrMax.size

    # take a view to have it 1d and not to mess with the plot
    # no copy is involved
    Mmax1d = Mmax.view()
    Mmax1d.shape = Mmax.size

    t=Table([Mmax1d, massRank, boxID, fileNrMax, subNrMax],
            names=['M200', 'MassRank', 'boxID', 'fileNr', 'subNr'])

    t['M200'].unit = u.dex(u.Msun)
    t['MassRank'].description = ('Rank of the halo in mass sorting'
                        f" ([0..{args.n-1}], 0 is the most massive)")
    t['boxID'].description = ('ID of the box used for bootstrap, '
                        f" [0..{args.Nboxes}]")

    t.meta = {'Description':'Hebb result table','Nboxes':args.Nboxes,
              'z_target':args.z_target, '-n':args.n,
              '--survey':args.survey, '-L':args.L, '-M': args.M}

    t.sort(['MassRank', 'boxID', 'fileNr', 'subNr'])
    t.write('hebb_halo_samples.txt', format='ascii.ecsv', overwrite=True)

    loggerH(f"TABLE: written 'hebb_halo_samples.txt'")




def pretty_plot(Mmax: npt.NDArray) -> None:
    """
    Pretty plot the distribution of the bootstrap results.

    Parameters
    ----------
    Mmax : np.ndarray of shape (NMassRanks, Nboxes)
        A 2D array of floats representing the maximum halo masses
        found in each bootstrap iteration. Units should be log10(M_sun).

    Returns
    -------
    None

    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    Mmaxmin = Mmax.min()
    Mmaxmax = Mmax.max()
    for i in range(Mmax.shape[0]):
        hist, bins = np.histogram(Mmax[i,:], 50, range=[Mmaxmin, Mmaxmax])
        ax.plot((bins[:-1]+bins[1:])/2, hist, label=f"{i}")
    ax.set_ylabel('N halos')
    ax.set_xlabel(r'$\log (M_{200}/M_\odot)$')
    ax.legend(loc='best')
    plt.show()
