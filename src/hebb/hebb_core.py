#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import numpy.typing as npt
import h5py
from pathlib import PurePath, Path
import astropy.units as u
import astropy.cosmology.units as cu
u.add_enabled_units(cu)
from astropy.units import Quantity
import numpy as np
# from numba import njit
from scipy import spatial
from time import perf_counter
from .logger_mine import loggerH, loggerN



def comovSideLenght(area: Quantity, zmin: float , zmax: float) -> Quantity:
    """
    Compute the side lenght of a cube of a volume equal to the estimated volume
    of a survey.

    Parameters
    ----------
    area: astropy.units.Quantity
        Area of the survey, must be an astropy unit quantity, e.g.
        300*(u.arcmin**2).

    zmin: float
        Min redshift of the survey.

    zmax: float
        Max redshift of the survey.

    Returns
    -------
    L: astropy.units.Quantity
        Side lenght of a box in cMpc.
    """

    from astropy.cosmology import Planck15
    Omega     = area.to(u.steradian).value # get rid of unit
    d2        = Planck15.comoving_distance(zmin)
    d3        = Planck15.comoving_distance(zmax)
    V         = Omega/3 * (d3**3 - d2**3)
    newL      = np.cbrt(V.to(u.Mpc**3))

    return newL


def loadCatalogue(z_target: float,
                  path_data: PurePath | Path | str,
                  *,
                  M: float | None = None,
                  force_light: bool = False) -> tuple[npt.NDArray[np.float32],
                                                      npt.NDArray[np.float32],
                                                      npt.NDArray[np.float32],
                                                      float]:
    """
    Load the Uchuu catalogue, possibly using a mass cut.

    Parameters
    ----------
    z_target: float
        Redshift of the target.

    path_data: string or pathlib object
        Path of the Uchuu simulation catalogue.

    M: float, optional
        Database mass cut in Msun, greatly speed up the database loading and search but
        you can incur into empty selection. Default None (load all the
        database).

    force_light: bool, optional
        Force the reading of the light catalogue first. Default False.



    Returns
    -------
    coords: np.ndarray of shape(NSubHaloes, 3)
        Subhalo coordinates from catalogue in cMpc.

    centers: np.ndarray of shape(Nboxes, 3)
        Center of boxes used for bootstrap.

    mass: np.ndarray of shape(NSubHaloes,)
        Subhalo masses in Msun.

    fileNr: np.ndarray of shape(NSubHaloes,)
        Merger tree file number of the subhaloes, it is sampled and given in
        output for tracing a subhalo with the merger tree.

    subNr: np.ndarray of shape(NSubHaloes,)
        Subhalo number in a specific merger tree file, it is sampled and given in
        output for tracing a subhalo with the merger tree.

    BoxSize: float
        Size of the box of the numerical simulation in cMpc.


    """

    from .uchuu_snaps_z import uchuu_snap_list
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
            if (np.log10(M) < m200_bins_edges[0] and not
                np.allclose(np.log10(M) , m200_bins_edges[0], atol=1e-3)):
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
    loggerN(f"CATALOGUE: read {mass.size} haloes")
    loggerN(f"CATALOGUE: log(M200/Msun): min={mass[0]:.2f},"
          f" max={mass[-1]:.2f}")

    return coords, mass, fileNr, subNr, BoxSize

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

def bootstrap_kdtree_single(BoxSize:    float,
                            newL:       float,
                            coords:     npt.NDArray[np.float32],
                            centers:    npt.NDArray[np.float32],
                            mass:       npt.NDArray[np.float32],
                            fileNr:     npt.NDArray[np.float32],
                            subNr:      npt.NDArray[np.float32],
                            *,
                            NMassRank: int = 1,
                            leafsize:  int = 128,
                            workers:   int = -1) -> tuple[npt.NDArray[np.float32],
                                                          npt.NDArray[np.float32],
                                                          npt.NDArray[np.float32]]:
    """
    Do a search using a single KDTree and perform a halo search within N boxes. Better
    than the double one since the centers are uniformly seeded in the volume,
    and it can be parallelized

    Parameters
    ----------
    BoxSize: float
        Size of the box of the numerical simulation in cMpc.

    newL: float
        Half size of the box to bootstrap in cMpc. This is used as the maximum
        (Chebyshev) distance to find the subhaloes during a single search

    coords: np.ndarray of shape(NSubHaloes, 3)
        Subhalo coordinates from catalogue in cMpc.

    centers: np.ndarray of shape(Nboxes, 3)
        Center of boxes used for bootstrap.

    mass: np.ndarray of shape(NSubHaloes,)
        Subhalo masses in Msun.

    fileNr: np.ndarray of shape(NSubHaloes,)
        Merger tree file number of the subhaloes, it is sampled and given in
        output for tracing a subhalo with the merger tree.

    subNr: np.ndarray of shape(NSubHaloes,)
        Subhalo number in a specific merger tree file, it is sampled and given in
        output for tracing a subhalo with the merger tree.

    NMassRank: int, optional
        Max rank in mass to search. The bootstrap will track the NMassRank most
        massive galaxies in each box. Default 1

    leafsize: int, optional
        Size of each leaf for the KDTree, increasing it speeds up the research.
        Default: 128

    workers: int, optional
        Number of jobs to schedule for parallel processing during the tree
        queries. If -1 is given all processors are used. Default: -1


    Returns
    -------

    Mmax: np.ndarray of shape(NMassRank, Nboxes)
        2D array of log10(M200/Msun) masses for all the mass ranks

    fileNrMax: np.ndarray of shape(NMassRank, Nboxes)
        2D array of merger tree file IDs for each halo for traceback

    subNrMax: np.ndarray of shape(NMassRank, Nboxes)
        2D array of subhalo IDs in each of the merger tree files for each halo for traceback

    """

    # Number of boxes for the bootstrap
    Nboxes = centers.shape[0]

    Mmax = np.zeros((NMassRank,Nboxes), dtype=mass.dtype)
    fileNrMax = np.zeros((NMassRank,Nboxes), dtype=fileNr.dtype)
    subNrMax = np.zeros((NMassRank,Nboxes), dtype=subNr.dtype)

    # unfortunately as for scipy 1.17.1, only float64 internal operations are
    # supported, the documentation says that everything will be copied
    # internally to float64, even with copy_data=False.
    # There is a request for a float32 support https://github.com/scipy/scipy/issues/18467
    # hopefully it will get implemented. This is the only bottleneck in the
    # code when you have over 100 millions haloes
    t0 = perf_counter()
    loggerH(f"KDtree: Building KDtree for fast search")
    tree = spatial.KDTree(coords, boxsize=BoxSize, leafsize=leafsize,
                          copy_data=False)
    loggerN(f"KDtree: building time {perf_counter()-t0:.1f} s")

    loggerH(f"SEARCH: starting {Nboxes} iterations...")
    t0 = perf_counter()

    dtype=np.min_scalar_type(coords.shape[0])

    # keeping the loop in pure python is ok, the bottleneck is the query
    for j in range(Nboxes):
        ind = np.array(tree.query_ball_point(centers[j,:], newL, p=np.inf,
                                             workers=workers, return_sorted = False),
                       dtype=dtype)

        if ind.size < NMassRank:
            raise ValueError('Hitting a region with less than '
                             f"{NMassRank} subhalos, in case you used -M try"
                             ' to lower (or omit) the mass cut')

        mass_tmp = mass[ind]
        fileNr_tmp = fileNr[ind]
        subNr_tmp = subNr[ind]

        index = np.argsort(mass_tmp)[-NMassRank:][::-1]
        Mmax[:,j] = mass_tmp[index]
        fileNrMax[:,j] = fileNr_tmp[index]
        subNrMax[:,j] = subNr_tmp[index]

    loggerN(f"SEARCH: done, time {perf_counter()-t0:.1f} s")

    return Mmax, fileNrMax, subNrMax



def hebb(z_target: float,
         Nboxes: int,
         path_data: PurePath | Path | str,
         *,
         NMassRank: int = 1,
         survey: tuple[float,float,float] | None = None,
         L: float | None = None,
         M: float | None = None,
         method: str = 'bootstrap',
         leafsize: int = 128,
         kdtreeWorkers: int = -1,
         force_light: bool = False) -> tuple[npt.NDArray[np.float32],
                                             npt.NDArray[np.float32],
                                             npt.NDArray[np.float32]]:
    """
    Compute the N most massive dark matter haloes that you can find in a given
    survey with [z_min, z_max] and field-of-view, or by directly passing a box
    size, by performing a non-parametric block bootstrap over the Uchuu (2/h cGpc)^3 run.

    Parameters
    ----------
    z_target: float
        Redshift of the target.

    Nboxes: int
        Number of boxes (iterations) for the block bootstrap. Needs to be a
        positive integer.

    path_data: string or pathlib object
        Path of the Uchuu simulation catalogue.

    NMassRank: int, optional
        Max rank in mass to search. The bootstrap will track the NMassRank most
        massive galaxies in each box. Default: 1 (only the most massive halo is
        tracked).

    survey: tuple of floats, (z_min, z_maz, fov), optional
        Survey min z, max z, FOV in arcmin^2, used to estimate the box size for
        the block bootstrap. Mutually exclusive with L. Default: None.

    L: float, optional
        Size of the box for block bootstrap in cMpc. Mutually exclusive with
        `survey`. Default: None.

    M: float, optional
        Database mass cut in Msun, greatly speed up the database loading and search but
        you can incur into empty selection. Default: None (load all the
        database).

    leafsize: int, optional
        Size of each leaf for the KDTree, increasing it speeds up the research.
        Default: 128.

    kdtreeWorkers: int = -1,
        Number of jobs to schedule for parallel processing during the tree
        queries. If -1 is given all processors are used. Default: -1

    force_light: bool, optional
        Force the reading of the light catalogue first. Default False.

    Returns
    -------

    Mmax: np.ndarray of shape(NMassRank, Nboxes)
        2D array of log10(M200/Msun) masses for all the mass ranks

    fileNrMax: np.ndarray of shape(NMassRank, Nboxes)
        2D array of merger tree file IDs for each halo for traceback

    subNrMax: np.ndarray of shape(NMassRank, Nboxes)
        2D array of subhalo IDs in each of the merger tree files for each halo for traceback

    """



    if survey is None and L is None:
        raise ValueError('L and survey cannot be both None')

    coords, mass, fileNr, subNr, BoxSize = loadCatalogue(z_target, path_data,
                                                         force_light = force_light,
                                                         M = M)


    # L has a higher priority over survey
    if L is None:
        # this code is here just for mypy, it's alreay checked at the beginning
        if survey is None:
            raise ValueError('L and survey cannot be both None')

        # for surveys, (zmin,zmax,fov)
        area = survey[2]*(u.arcmin**2)
        newL = comovSideLenght(area,survey[0], survey[1]).value/2 # in cMpc
        loggerH(f"BOX: Estimated volume={8*newL**3:.3e} cMpc^3,  Box size "
                f"L={newL*2:.3f} cMpc")
    else:
        newL = L/2.
        loggerH(f"BOX: Volume={8*newL**3:.3e} cMpc^3,  Box size L={newL*2:.3f} cMpc")



    rng = np.random.default_rng()
    if method == 'montecarlo':
        # shoot (Nboxes,3) random numbers between 0 and BoxSize
        centers=rng.random(size=(Nboxes,3), dtype=np.float32)
        centers*=BoxSize
    elif method == 'bootstrap':
        # tiling, no overlapping
        x = np.arange(0., BoxSize, newL, dtype=np.float32)
        xc = (x[1:]+x[:-1])/2.

        # move a bit the grid so we randomly get everything
        off=rng.random(size=(1), dtype=np.float32)
        off*=BoxSize - x[-1]
        xc += off

        XC, YC, ZC = np.meshgrid(xc, xc, xc, copy=True, sparse=False)
        XC.shape = XC.size
        YC.shape = YC.size
        ZC.shape = ZC.size
        centers = np.zeros((XC.size, 3), dtype=XC.dtype)
        centers[:,0] = XC
        centers[:,1] = YC
        centers[:,2] = ZC


    # Mmax, fileNrMax, subNrMax = bootstrap_brute_force(Nboxes, BoxSize, newL, coords,
                                          # centers, mass, fileNr, subNr)
    # Mmax, fileNrMax, subNrMax = bootstrap_kdtree_double(Nboxes, BoxSize, newL, coords,
                                          # centers, mass, fileNr, subNr)
    Mmax, fileNrMax, subNrMax = bootstrap_kdtree_single(BoxSize, newL, coords,
                                                        centers, mass, fileNr, subNr,
                                                        leafsize=leafsize,
                                                        NMassRank=NMassRank,
                                                        workers=kdtreeWorkers)

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
