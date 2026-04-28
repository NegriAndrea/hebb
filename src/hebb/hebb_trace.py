#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from astropy.table import Table, vstack
from .logger_mine import loggerH, loggerN

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

def chopArray(arraySize, Nchunks):
    """
    Code taken from array_split, except that it takes the size of an array, and
    returns the offsets and dimensions of the subarrays, divided in Nchunks of
    approximately equal size.

    arraySize: size of the array to chop
    Nchunks : number of chunks to dive the array in

    Works on 1D array only.

    Returns:

    offsets : (Nchunks,) integer array
    sizes   : (Nchunks,) integer array

    """
    import numpy as np

    # N is a scalar
    Nsections = int(Nchunks)

    if Nsections <= 0:
        raise ValueError('N must be larger than 0.')
    Neach_section, extras = divmod(arraySize, Nsections)
    sizes = (extras * [Neach_section+1] +
                     (Nsections-extras) * [Neach_section])

    offsets = np.array([0]+sizes[:-1], dtype=np.intp).cumsum()
    sizes = np.asarray(sizes, dtype=np.intp)

    assert offsets.shape == (Nchunks,)
    assert sizes.shape == (Nchunks,)
    assert offsets[-1]+sizes[-1] == arraySize

    return offsets, sizes

def hebb_trace(tableName, mergerTreePath, targetZ, v, serial):

    if serial:
        # do it serially
        t = Table.read(tableName, format='ascii.ecsv')
        return (hebb_trace_single(tableName, mergerTreePath, targetZ, v, 0,
                len(t)),0)
    else:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        t = Table.read(tableName, format='ascii.ecsv')
        t.sort(['fileNr'])

        offset, sizes = chopArray(len(t), size)

        i_begin = offset[rank]
        i_end = sizes[rank]+i_begin

        local_t = hebb_trace_single(tableName,
                mergerTreePath, targetZ, v, i_begin, i_end)

        all_tables = comm.gather(local_t, root=0)

        if rank == 0:
            final_table = vstack(all_tables)
            return final_table, rank
        else:
            return None, rank




def hebb_trace_single(tableName, mergerTreePath, targetZ, v, i_begin, i_end,
                      Ntrack = 10):
    """
    Trace back galaxies found with hebb.

    """
    import numpy as np
    from pathlib import PurePath
    import astropy.units as u
    import astropy.cosmology.units as cu
    u.add_enabled_units(cu)
    import chydrotree

    from .uchuu_snaps_z import uchuu_snap_list

    snapNr_list, z = uchuu_snap_list()
    targetSnapNr = snapNr_list[np.abs(z - targetZ).argmin()]

    t = Table.read(tableName, format='ascii.ecsv')
    t=t[i_begin:i_end]

    t.sort(['fileNr', 'subNr'])

    UfileNr, off, size = unique_ordered(t['fileNr'])

    if v>1:
        tmp = Table([UfileNr, size], names=['MergTreeID', 'Nhalos'])
        print(tmp)
        del tmp

    fullT = []
    fileNr_original = []
    subNr_original = []
    massMinorMerger = []
    massSingleProgs = []
    NMinorProgs = []

    tmp = np.full(Ntrack, -np.inf, dtype=np.float32)

    # I loop over the merger tree files (MT)
    for jMT, o, s in zip(UfileNr, off, size):

        # table containing only files related to 1 merger tree file
        subt = t[o:o+s]

        fileMT = PurePath(mergerTreePath) / f'mergertree_{jMT}.h5'

        loggerN(f'Read {fileMT}')
        # in case I want the position, use newFields=['Pos']
        tree = chydrotree.forestCT(fileMT, newFields=[])

        for sNr in subt['subNr']:
            loggerN(f"Doing {sNr}")

            try:
                # depends what I want to find, in this case the progenitor
                # do I need raiseError=True?
                tHist = tree.subHistory(sNr, targetSnapNr=targetSnapNr,
                                        t=True, raiseError=True)

                # filter to get the target I want
                tHist = tHist[tHist['SnapNr'] == targetSnapNr]
                fullT.append(tHist)
                fileNr_original.append(jMT)
                subNr_original.append(sNr)

                tall = tree.subAllProgs(sNr, t=True)
                tall = tall[tall['SnapNr'] == targetSnapNr]

                # filter out the main progenitor
                tall = tall[tall['SubhaloNr'] != tHist['SubhaloNr']]
                tall.sort('Mass')
                tall = tall[::-1]

                massMinorMerger.append(np.sum(tall['Mass'].to('Msun')).value)
                NMinorProgs.append(len(tall))

                # track the first Ntrack haloes
                tmp[:] = -np.inf
                tmp[:min(len(tall), tmp.size)] = (
                        tall['Mass'][:min(len(tall), tmp.size)])
                massSingleProgs.append(tmp)
            except chydrotree.CTNoProgenitorError:
                if v>1:
                    print('NoProgenitor')

    fullT = vstack(fullT)
    tmp = np.vstack(massSingleProgs)
    assert tmp.shape[1] == Ntrack

    fullT['massMinorProgs'] = np.log10(np.array(massMinorMerger))
    fullT['massMinorProgs'].unit = u.dex('Msun')
    fullT['NMinorProgs'] = np.array(NMinorProgs)

    for i in range(Ntrack):
        name = str(i+1)+'prog'
        fullT[name] = tmp[:,i]
        fullT[name].unit = u.dex('Msun')

    fullT['fileNr_original'] = np.array(fileNr_original,
            dtype=t['fileNr'].dtype)
    fullT['subNr_original'] = np.array(subNr_original,
            dtype=t['subNr'].dtype)

    return fullT
