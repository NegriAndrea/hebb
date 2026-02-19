#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def unique_ordered(x):
    """
    Same output of np.unique when return_counts=True,
    return_index=True, but 1000 faster since x is assumed sorted.
    It works on 1d arrays, untested for other shapes.
    """
    import numpy as np

    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError('x must be a 1D array-like')

    # compare the elements and find where the elements change. This
    # is better than using np.diff, since the comparison always
    # work, while for some arrays (like recarrays) the difference
    # may be not defined
    flt =  np.flatnonzero(x[1:] != x[:-1])

    # by construction, we to add 1 to flt
    off = np.concatenate([[0], flt+1])
    dim = np.concatenate([np.diff(off), [x.size-off[-1]]])
    unique = x[off]

    return unique, off, dim

def hebb_trace(tableName, mergerTreePath, targetZ):
    """
    Trace back galaxies found with hebb.

    """
    import numpy as np
    from pathlib import PurePath
    from astropy.table import Table, vstack
    import astropy.units as u
    import astropy.cosmology.units as cu
    u.add_enabled_units(cu)
    import helperspy as hpy

    from .uchuu_snaps_z import uchuu_snap_list

    snapNr_list, z = uchuu_snap_list()
    targetSnapNr = snapNr_list[np.abs(z - z_target).argmin()]

    t = Table.read(tableName, format='ascii.ecsv')

    t.sort(['fileNr', 'subNr'])

    UfileNr, off, size = unique_ordered(t['fileNr'])

    fullT = []
    fileNr_original = []
    subNr_original = []

    # I loop over the merger tree files (MT)
    for jMT, o, s in zip(UfileNr, off, size):

        # table containing only files related to 1 merger tree file
        subt = t[o:o+s]

        fileMT = PurePath(mergerTreePath) / f'mergerTree_{jMT}.hdf5'

        # in case I want the position, use newFields=['Pos']
        tree = hpy.forestCT(fileMT, newFields=[])

        for sNr in subt['subNr']:

            try:
                # depends what I want to find, in this case the progenitor
                # do I need raiseError=True?
                tHist = tree.subHistory(sNr, targetSnapNr=targetSnapNr,
                                        t=True, raiseError=True)

                # filter to get the target I want
                fullT.append(tHist[tHist['snapNr'] == targetSnapNr])
                fileNr_original.append(jMT)
                subNr_original.append(sNr)
            except hpy.NoProgenitorError:
                pass

    fullT = vstack(fullT)
    fullT['fileNr_original'] = np.concatenate(fileNr_original)
    fullT['subNr_original'] = np.concatenate(subNr_original)

    return fullT
