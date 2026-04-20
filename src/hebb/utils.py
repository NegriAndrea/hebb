#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import numpy.typing as npt
from .logger_mine import loggerH, loggerN
from astropy.table import Table

def pretty_print(Mmax: npt.NDArray,
                 *,
                 write: bool = False,
                 tvar: npt.NDArray | None = None) -> None:
    """
    Pretty print a summary of the results and quantile table.

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
    NMassRanks = Mmax.shape[0]
    import astropy.units as u
    quantiles = np.array([0.16, 0.5, 0.86,
                          0.05, 0.1, 0.25, .75, 0.9, 0.95])
    quants = np.zeros((Mmax.shape[0], quantiles.size))

    for i in range(Mmax.shape[0]):
        quants[i,:] = np.quantile(10.**Mmax[i,:].astype(np.float64),
                                            quantiles)
        M16, M50, M86 = quants[i,:3]
        # print(f"{i} log10(M200) = {np.log10(M50):.2f} _-{(M50-M16)/M50/np.log(10):.2f} ^+{(M86-M50)/M50/np.log(10):.2f}")

    # build a quantile table and use astropy Table to do a pretty print
    ind = [3,4,0,5,1,6,2,7,8]
    quantslog = np.log10(quants[:,ind])
    names=[str(int(q)) for q in quantiles[ind]*100]

    t2=Table(quantslog, names=names)

    # add MassRank column
    t2.add_column(np.arange(quantslog.shape[0], dtype=np.uint8), name='MassRank',
                  index=0)

    # print only 2 decimals for floating point numbers
    for cname in t2.colnames:
        if t2[cname].info.dtype in ['<f4', '<f8']:
            t2[cname].info.format = '6.2f'




    t3 = t2[['MassRank']]
    t3['log10(M200)'] = np.log10(quants[:,1])

    t3['+'] = (quants[:,1]-quants[:,0])/quants[:,1]/np.log(10)
    t3['-'] = (quants[:,2]-quants[:,1])/quants[:,1]/np.log(10)

    # print only 2 decimals for floating point numbers
    for cname in t3.colnames:
        if t3[cname].info.dtype in ['<f4', '<f8']:
            t3[cname].info.format = '6.2f'
    t3['log10(M200)'].units = u.dex('Msun')

    print('')
    print('')
    print('RESULTS: log median for every mass rank with ± 16th and 68th percentiles')
    print('         in latex is $log(M200)^{+}_{-}$')
    print(t3)

    print('')
    print('PERCENTILES TABLE')
    print(t2)

    if tvar is not None:
        if NMassRanks != tvar.shape[0]:
            raise ValueError('tvar does not contain the same number of mass'
                             ' ranks as the mass array')
        print('')
        print('PERCENTILES VARIANCE')
        print('i.e. how much we trust the values of the percentiles table')
        tvar = Table(tvar, copy=False, names = ('16', '50', '86'))
        tvar.add_column(np.arange(NMassRanks, dtype=np.uint8), name='MassRank',
                     index=0)
        for cname in tvar.colnames:
            if tvar[cname].info.dtype in ['<f4', '<f8']:
                tvar[cname].info.format = '6.1e'
        print(tvar)


    if write:
        t2.write('hebb_percentiles.txt', format='ascii.ecsv', overwrite=True)
        loggerH(f"TABLE: written 'hebb_quantiles.txt'")



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
                        f" [0..{Mmax.shape[1]}]")

    t.meta = {'Description':'Hebb result table','Nboxes':boxID.max(),
              'z_target':args.z_target, '-n':args.n,
              '--survey':args.survey, '-L':args.L, '-M': args.M}

    t.sort(['MassRank', 'boxID', 'fileNr', 'subNr'])
    t.write('hebb_halo_samples.txt', format='ascii.ecsv', overwrite=True)

    loggerN(f"TABLE: written 'hebb_halo_samples.txt'")




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
        hist, bins = np.histogram(Mmax[i,:], 50, range=(Mmaxmin, Mmaxmax))
        ax.plot((bins[:-1]+bins[1:])/2, hist, label=f"{i}")
    ax.set_ylabel('N halos')
    ax.set_xlabel(r'$\log (M_{200}/M_\odot)$')
    ax.legend(loc='best')
    plt.show()
