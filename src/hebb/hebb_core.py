#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import h5py
from pathlib import PurePath, Path
import astropy.units as u
import astropy.cosmology.units as cu
u.add_enabled_units(cu)
import numpy as np
# from numba import njit
from scipy import spatial
import logging
from time import perf_counter

logger = logging.getLogger(__name__)

def loggerN(msg):
    logger.info(msg)

def loggerH(msg):
    logger.info('')
    logger.info(msg)

def comov_volume(area, z1, z2):
    """
    Area must be with an astropy unit, e.g. 300*(u.arcmin**2)
    Returns size lenght of a box in cMpc

    """
    from astropy.cosmology import Planck15
    Omega     = area.to(u.steradian).value # get rid of unit
    d2        = Planck15.comoving_distance(z1)
    d3        = Planck15.comoving_distance(z2)
    V         = Omega/3 * (d3**3 - d2**3)
    newL      = np.cbrt(V.to(u.Mpc**3))

    return newL

# @njit
def dL(coord1, coord2_in, boxsize):
    """
    Computes the distance in a 3D periodic box between a set of coordinates and
    a point, or another, broadcastable set of coordinates.

    coord1: array of dimensions (n,3)

    coord2: array of dimensions (1,3) or (3,) in the case of a single point, or
                array of size (n,3) in case of multiple points

    boxsize : scalar, the size of the box, each dimension has the same size

    Returns: the distance squared


    usage:

    import numpy as np

    box = 0.5

    # generate fake data
    rng = np.random.default_rng(12345)
    c1 = rng.random((1000,3))

    distancePeriodic(c1, np.array([0.,0.,0.], box))


    c2 = rng.random((1000,3))
    distancePeriodic(c1, c2, box))

    """

    message = 'coord1 must be a two dimensional array with shape[1]==3'
    coord2=np.atleast_2d(coord2_in)

    if coord1.ndim != 2:
        raise ValueError(message)

    if coord1.shape[1] != 3:
        raise ValueError(message)

    if coord2.ndim == 1:
        if coord2.size != 3:
            raise ValueError('I need 3 points as np.array([x,y,z])')
    elif coord2.ndim == 2:
        if coord2.shape[1] != 3:
            raise ValueError('I need 3 points as np.array([x,y,z])')
    else:
        raise ValueError('I need 3 points as np.array([x,y,z])')

    dx = coord1[:,0] - coord2[:,0]
    dy = coord1[:,1] - coord2[:,1]
    dz = coord1[:,2] - coord2[:,2]

    dx[(dx > boxsize*0.5)] -= boxsize
    dx[(dx < - boxsize*0.5)] += boxsize
    dy[(dy > boxsize*0.5)] -= boxsize
    dy[(dy < - boxsize*0.5)] += boxsize
    dz[(dz > boxsize*0.5)] -= boxsize
    dz[(dz < - boxsize*0.5)] += boxsize

    return dx, dy, dz

# @njit(parallel=True)
def bootstrap_brute_force(Nboxes, BoxSize, newL, coords, centers, mass, fileNr, subNr):

    Mmax = np.zeros(Nboxes)
    fileNrMax = np.zeros(Nboxes, dtype=fileNr.dtype)
    subNrMax = np.zeros(Nboxes, dtype=subNr.dtype)

    for j in range(Nboxes):

        dx, dy, dz = dL(coords, centers[j,:], BoxSize)
        mask = (np.abs(dx) < newL) & (np.abs(dy) < newL) & (np.abs(dz) < newL)

        if np.count_nonzero(mask) == 0:
            raise ValueError('Hitting an empty region, in case you used -M try'
                             ' to lower (or omit) the mass cut')
        mass_tmp = mass[mask]
        fileNr_tmp = fileNr[mask]
        subNr_tmp = subNr[mask]

        index = np.argmax(mass_tmp)

        Mmax[j] = mass_tmp[index]
        fileNrMax[j] = fileNr_tmp[index]
        subNrMax[j] = subNr_tmp[index]

    return Mmax, fileNrMax, subNrMax

def bootstrap_kdtree_double(Nboxes, BoxSize, newL, coords, centers, mass, fileNr, subNr):
    """
    Do a search using a double KDTree and perform a block bootstrap.

    """

    Mmax = np.zeros(Nboxes)
    fileNrMax = np.zeros(Nboxes, dtype=fileNr.dtype)
    subNrMax = np.zeros(Nboxes, dtype=subNr.dtype)

    tree_data = spatial.KDTree(coords, boxsize=BoxSize, leafsize=32)
    tree_centres = spatial.KDTree(centers, boxsize=BoxSize, leafsize=32)
    indexes = tree_centres.query_ball_tree(tree_data, newL, p=np.inf)

    assert len(indexes) == Nboxes
    for j in range(Nboxes):

        if len(indexes[j]) == 0:
            raise ValueError('Hitting an empty region, in case you used -M try'
                             ' to lower (or omit) the mass cut')
        ind = np.array(indexes[j],
                             dtype=np.min_scalar_type(coords.shape[0]))
        mass_tmp = mass[ind]
        fileNr_tmp = fileNr[ind]
        subNr_tmp = subNr[ind]

        index = np.argmax(mass_tmp)
        Mmax[j] = mass_tmp[index]
        fileNrMax[j] = fileNr_tmp[index]
        subNrMax[j] = subNr_tmp[index]

    return Mmax, fileNrMax, subNrMax

def bootstrap_kdtree_single(Nboxes, BoxSize, newL, coords, centers, mass,
                            fileNr, subNr, leafsize=128, nn=1):
    """
    Do a search using a single KDTree and perform a block bootstrap. Better
    than the double one since the centers are uniformly seeded in the volume,
    and it can be parallelized

    """

    Mmax = np.zeros((nn,Nboxes), dtype=mass.dtype)
    fileNrMax = np.zeros((nn,Nboxes), dtype=fileNr.dtype)
    subNrMax = np.zeros((nn,Nboxes), dtype=subNr.dtype)

    t0 = perf_counter()
    loggerH(f"KDtree: Building KDtree for fast search")
    tree = spatial.KDTree(coords, boxsize=BoxSize, leafsize=leafsize)
    loggerN(f"KDtree: building time {perf_counter()-t0:.1f} s")

    loggerH(f"BLOCK BOOTSTRAP: startig...")
    t0 = perf_counter()

    dtype=np.min_scalar_type(coords.shape[0])

    for j in range(Nboxes):
        ind = np.array(tree.query_ball_point(centers[j,:], newL, p=np.inf, workers=-1), dtype=dtype)

        if ind.size < nn:
            raise ValueError('Hitting a region with less than '
                             f"{nn} subhalos, in case you used -M try"
                             ' to lower (or omit) the mass cut')

        mass_tmp = mass[ind]
        fileNr_tmp = fileNr[ind]
        subNr_tmp = subNr[ind]

        index = np.argmax(mass_tmp)
        index = np.argsort(mass_tmp)[-nn:][::-1]
        Mmax[:,j] = mass_tmp[index]
        fileNrMax[:,j] = fileNr_tmp[index]
        subNrMax[:,j] = subNr_tmp[index]

    loggerN(f"BLOCK BOOTSTRAP: done, time {perf_counter()-t0:.1f} s")

    return Mmax, fileNrMax, subNrMax



def hebb(z_target, Nboxes, path_data, *, survey=None, L=None, M=None,
         leafsize=128, force_light=False, nn = 1):
    from .uchuu_snaps_z import uchuu_snap_list

    if survey is None and L is None:
        raise ValueError('L and survey cannot be both None')

    snapNr_list, z = uchuu_snap_list()
    snapNr = snapNr_list[np.abs(z - z_target).argmin()]

    # Boxsize in cMpc
    BoxSize = 2000/0.6774

    catFileName = Path(path_data)/f'catalogue_uchuu.hdf5'
    if not catFileName.is_file() or force_light:
        catFileName = Path(path_data)/f'catalogue_uchuu_light.hdf5'

    if not catFileName.is_file():
        raise IOError('I cannot locate the catalogue file, I have tried'
                      f"{Path(path_data)/f'catalogue_uchuu.hdf5'}"
                      f" and {Path(path_data)/f'catalogue_uchuu_light.hdf5'}"
                      f" visit https://github.com/NegriAndrea/hebb")

    loggerH(f"CATALOGUE: Reading {catFileName}")
    t0 = perf_counter()

    with h5py.File(catFileName, 'r') as ff:

        # the catalogue is sorted in M200, with histograms to load only the
        # haloes above a certain mass without having to read the full dataset
        # if not needed
        if M is None:
            offset = 0
        else:
            m200_indexes = ff[f'S-{snapNr}/M200_indexes'][()]
            m200_bins_edges = ff[f'S-{snapNr}/M200c_bins_edges'][()]
            if np.log10(M) < m200_bins_edges[0]:
                raise ValueError(f"The requested mass {M=:.2e} Msun is too low for "
                                 f"z={z_target}, database min mass at this z "
                                 f"is {10.**float(m200_bins_edges[0]):.2e} Msun")
            tmp = np.searchsorted(m200_bins_edges, np.log10(M))
            offset = int(m200_indexes[max(tmp-1,0)])

        coords = ff[f'S-{snapNr}/Coordinates'][offset:]
        mass= ff[f'S-{snapNr}/M200c'][offset:]
        fileNr = ff[f'S-{snapNr}/fileNr'][offset:]
        subNr = ff[f'/S-{snapNr}/SubNr'][offset:]
        bin_size = BoxSize/65_535

    coords=coords.astype(np.float32)*bin_size

    loggerN(f"CATALOGUE: reading time {perf_counter()-t0:.1f} s")
    loggerN(f"CATALOGUE: log(M200/Msun): min={mass[0]:.2f},"
          f" max={mass[-1]:.2f}")

    # L has a higher priority over survey
    if L is None:
        # for surveys
        area = survey['fov']*(u.arcmin**2)
        newL = comov_volume(area,survey['zmin'], survey['zmax']).value/2 # in cMpc
        loggerH(f"BOX: Estimated volume={8*newL**3:.3e} cMpc^3,  Box size "
                f"L={newL*2:.3f} cMpc")
    else:
        newL = L/2.
        loggerH(f"BOX: Volume={8*newL**3:.3e} cMpc^3,  Box size L={newL*2:.3f} cMpc")



    # shoot (Nboxes,3) random numbers between 0 and BoxSize
    rng = np.random.default_rng()
    centers=rng.random(size=(Nboxes,3), dtype=np.float32)
    centers*=BoxSize


    # Mmax, fileNrMax, subNrMax = bootstrap_brute_force(Nboxes, BoxSize, newL, coords,
                                          # centers, mass, fileNr, subNr)
    # Mmax, fileNrMax, subNrMax = bootstrap_kdtree_double(Nboxes, BoxSize, newL, coords,
                                          # centers, mass, fileNr, subNr)
    Mmax, fileNrMax, subNrMax = bootstrap_kdtree_single(Nboxes, BoxSize, newL, coords,
                                                        centers, mass, fileNr, subNr,
                                                        leafsize=leafsize,
                                                        nn=nn)

    return Mmax, fileNrMax, subNrMax

def hebb_estimate(z_target, path_data, L, M=None, v=0):
    from .uchuu_snaps_z import uchuu_snap_list

    snapNr_list, z = uchuu_snap_list()
    snapNr = snapNr_list[np.abs(z - z_target).argmin()]

    # Boxsize in cMpc
    BoxSize = 2000/0.6774

    catFileName = Path(path_data)/f'catalogue_uchuu.hdf5'
    if not catFileName.is_file():
        catFileName = Path(path_data)/f'catalogue_uchuu_light.hdf5'

    if not catFileName.is_file():
        raise IOError('I cannot locate the catalogue file, I have tried'
                      f"{Path(path_data)/f'catalogue_uchuu.hdf5'}"
                      f" and {Path(path_data)/f'catalogue_uchuu_light.hdf5'}")

    loggerN(f"Reading {catFileName}")

    with h5py.File(catFileName, 'r') as ff:

        # the catalogue is sorted in M200, with histograms to load only the
        # haloes above a certain mass without having to read the full dataset
        # if not needed
        if M is None:
            offset = 0
        else:
            m200_indexes = ff[f'S-{snapNr}/M200_indexes'][()]
            m200_bins_edges = ff[f'S-{snapNr}/M200c_bins_edges'][()]
            if np.log10(M) < m200_bins_edges[0]:
                raise ValueError(f"The requested mass {M=:.2e} Msun is too low for "
                                 f"z={z_target}, database min mass at this z "
                                 f"is {10.**float(m200_bins_edges[0]):.2e} Msun")
            tmp = np.searchsorted(m200_bins_edges, np.log10(M))
            offset = int(m200_indexes[max(tmp-1,0)])

        coords = ff[f'S-{snapNr}/Coordinates'][offset:]
        mass= ff[f'S-{snapNr}/M200c'][offset:]
        fileNr = ff[f'S-{snapNr}/fileNr'][offset:]
        subNr = ff[f'/S-{snapNr}/SubNr'][offset:]
        bin_size = BoxSize/65_535

    coords=coords.astype(np.float32)*bin_size

    loggerN(f"READ log(M200/Msun): min={mass[0]:.2f},"
          f" max={mass[-1]:.2f}")


    for n in range(32):
        nbins = 2**n
        H, bins = np.histogramdd(coords, bins=nbins)
        stop_cond = np.any(H==0)
        print(f"L={BoxSize/nbins} {stop_cond}")
        if stop_cond:
            return
