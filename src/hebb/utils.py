#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import numpy.typing as npt
from .logger_mine import loggerH, loggerN
from astropy.table import Table
import astropy.units as u
import scipy.special

def listSnapshotTimes() -> Table:
    from .uchuu_snaps_z import uchuu_snap_list
    snapNr_list, z = uchuu_snap_list()

    t = Table([snapNr_list, z], names=['SnapNr', 'z'])
    return t


def buildPercentileTable(Mmax: npt.NDArray,
                         *,
                         units: bool = False) -> (Table, Table, Table):
    """
    Build a summary of the results as a quantile table.

    This function takes the extremal mass distribution and computes standard
    quantiles to provide a human-readable summary of the hebb analysis.

    Parameters
    ----------
    Mmax : np.ndarray of shape (NMassRanks, Nboxes)
        A 2D array of floats representing the maximum halo masses
        found in each bootstrap iteration. Units should be log10(M_sun).

    units : bool, optional
        Add units to the astropy table. Default False

    Returns
    -------
    tSigmas: astropy.table.Table
        Table containing the percentiles computed as gaussian equivalent
        sigmas.

    tRound: astropy.table.Table
        Table containing the percentiles computed with nice round numbers.

    tErrors: astropy.table.Table
        Table containing the 1 and 2 sigma errors.

    """
    NMassRanks = Mmax.shape[0]

    sigmaQuantiles = np.array([1,2,3,4,5])
    highSigmaQuant  = (0.5*(1+scipy.special.erf(sigmaQuantiles/np.sqrt(2))))
    lowSigmaQuant = (0.5*(1-scipy.special.erf(sigmaQuantiles/np.sqrt(2))))

    quantiles = np.array([0.05, 0.1, 0.25, 0.5, .75, 0.9, 0.95])

    full_quantiles = np.concatenate([quantiles, lowSigmaQuant, highSigmaQuant])
    quants = np.zeros((NMassRanks, full_quantiles.size))

    names=[str(int(q)) for q in quantiles*100]
    names_round = names.copy()
    names.extend([str(q)+'σ' for q in np.concatenate([-sigmaQuantiles,
                                                       sigmaQuantiles])])
    names = np.array(names)

    for i in range(Mmax.shape[0]):
        quants[i,:] = np.quantile(10.**Mmax[i,:].astype(np.float64),
                                            full_quantiles)

    ind = np.argsort(full_quantiles)
    quants = quants[:,ind]
    names = names[ind]
    t=Table(np.log10(quants), names=names)

    # add MassRank column
    t.add_column(np.arange(NMassRanks, dtype=np.uint8), name='MassRank',
                  index=0)

    # print only 2 decimals for floating point numbers
    for cname in t.colnames:
        if t[cname].info.dtype in ['<f4', '<f8']:
            t[cname].info.format = '6.2f'
        if units and cname != 'MassRank':
            t[cname].unit = u.dex('Msun')


    tErrors = t[['MassRank']]
    tErrors['log10(M200)'] = t['50']

    tErrors['-1σ'] = t['50'] - t['-1σ']
    tErrors['+1σ'] = t['1σ'] - t['50']
    tErrors['-2σ'] = t['50'] - t['-2σ']
    tErrors['+2σ'] = t['2σ'] - t['50']

    # print only 2 decimals for floating point numbers
    for cname in tErrors.colnames:
        if tErrors[cname].info.dtype in ['<f4', '<f8']:
            tErrors[cname].info.format = '6.2f'

    tmp = ['MassRank']
    tmp.extend(names_round)
    tRound = t[tmp]
    del tmp

    tSigmas = t[[str(q)+'σ' for q in np.concatenate([-sigmaQuantiles[::-1],
                                                       sigmaQuantiles])]]
    tSigmas.add_column(t['50'], name = '0σ', index=sigmaQuantiles.size)
    tSigmas.add_column(t['MassRank'], name = 'MassRank', index=0)

    return tSigmas, tRound, tErrors



def prettyPrint(tSigmas: Table,
                tRound: Table,
                tErrors: Table,
                *,
                write: bool = False,
                tvar: npt.NDArray | None = None) -> None:
    """
    Pretty print the percentiles tables.

    Parameters
    ----------
    tSigmas: astropy.table.Table
        Table containing the percentiles computed as gaussian equivalent
        sigmas.

    tRound: astropy.table.Table
        Table containing the percentiles computed with nice round numbers.

    tErrors: astropy.table.Table
        Table containing the 1 and 2 sigma errors.

    write: bool, optional
        Write the percentile table. Default False

    tvar: np.ndarray, optional
        Variance of the parameters estimated by the bayesian bootstrap. Default
        None

    Returns
    -------
    None

    """
    print('')
    print('')
    print('RESULTS: log median for every mass rank with ± 16th and 68th percentiles')
    print('         in latex is $log(M200)^{+}_{-}$')
    print(tErrors)

    print('')
    print('')
    print('PERCENTILES TABLE in round percentiles (dex(Msun))')
    print(tRound)


    print('')
    print('PERCENTILES TABLE in σ percentiles from the median (dex(Msun))')
    print(tSigmas)

    if tvar is not None:
        NMassRanks = tRound['MassRank'].max()+1
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
        for cname in tSigmas.colnames:
            if cname != 'MassRank':
                tSigmas[cname].unit = u.dex('Msun')
            tSigmas.rename_column(cname, cname.replace('σ', 'sig'))
        tSigmas.write('hebb_sigma_percentiles.txt', format='ascii.ecsv', overwrite=True)

        for cname in tRound.colnames:
            if cname != 'MassRank':
                tRound[cname].unit = u.dex('Msun')
            tRound.rename_column(cname, cname.replace('σ', 'sig'))

        tRound.write('hebb_round_percentiles.txt', format='ascii.ecsv', overwrite=True)
        loggerH(f"TABLE: written 'hebb_quantiles.txt'")



def save_table(Mmax: npt.NDArray,
               fileNrMax: npt.NDArray,
               subNrMax: npt.NDArray,
               args: argparse.Namespace,
               *,
               write: bool = False) -> Table:
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
    Mmax1d = np.reshape(Mmax1d, shape=(Mmax1d.size,), copy=False)

    t=Table([Mmax1d, massRank, boxID, fileNrMax, subNrMax],
            names=['M200', 'MassRank', 'boxID', 'fileNr', 'subNr'])

    t['M200'].unit = u.dex(u.Msun)
    t['MassRank'].description = ('Rank of the halo in mass sorting'
                        f" ([0..{args.n-1}], 0 is the most massive)")
    t['boxID'].description = ('ID of the box used for bootstrap, '
                        f" [0..{Mmax.shape[1]}]")

    t.meta = {'Description':'Hebb result table','Nboxes':boxID.max(),
              '-n':args.n,
              '--survey':args.survey, '-L':args.L}

    try:
        t.meta['z_target'] = args.z_target
    except AttributeError:
        pass

    try:
        t.meta['-M'] = args.M
    except AttributeError:
        pass

    t.sort(['MassRank', 'boxID', 'fileNr', 'subNr'])

    if write:
        t.write('hebb_halo_samples.txt', format='ascii.ecsv', overwrite=True)
        loggerN(f"TABLE: written 'hebb_halo_samples.txt'")

    return t




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
