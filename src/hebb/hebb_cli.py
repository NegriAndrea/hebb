#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .hebb_core import hebb, hebb_estimate, EmptyBoxError
from .hebb_trace import hebb_trace
from .utils import (save_table, pretty_plot, buildPercentileTable,
                    listSnapshotTimes, prettyPrint)
import argparse
import logging
from .bootstrap import bootstrap_driver
from astropy.table import vstack

def define_args() -> argparse.Namespace:
    """
    Define the command line interface and options of `hebb`.

    """

    description = (
            'Compute the N most massive dark matter haloes that you can find in a given '
            'survey with [z_min, z_max] and field-of-view with '
            'monte carlo, non-overlapping or non-parametric bayesian bootstrap over the Uchuu (2/h cGpc)^3 run.'
            )

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('z_target', type=float, help='Redshift of your target')

    group = parser.add_argument_group( "Processing mode (required, mutually exclusive)")
    group_container = group.add_mutually_exclusive_group(required=True)
    group_container.add_argument('--survey', nargs=3, type=float,
                                 help='Survey min z, max z, FOV in arcmin^2'
                                 ' (used to compute box volume)',
                                 metavar=("z_min", "z_max", "fov"))
    group_container.add_argument('-L', type=float, help='Size of the box in cMpc,'
                                 ' alternative to the boxsize computation from the FOV'
                                 ' and z-depth of the survey')

    parser.add_argument('-n', type=int, help='Track the N most massive haloes'
                        ' in each box (1 tracks only the most massive) [default: %(default)d]', default=1)
    parser.add_argument('--niter', type=int, help='Do a Monte Carlo selection'
                        ' with NITER overlapping boxes overlapping search',
                        default=0)
    parser.add_argument('--bb', type=int, help='Estimate the variance of the '
                        ' median with BB iterations bayesian non parametric '
                        ' bootstrap', default=0)
    parser.add_argument('-t', action='store_true', help='Create a table with'
                        ' the sampled haloes')
    parser.add_argument('--plot', action='store_true', help='Show a plot of the M200 distribution')
    parser.add_argument('-l', action='store_true', help='Show the list of '
                        ' snapshots redshifts and exit')

    parser.add_argument('-M', type=float, help='OPTIMIZATION: Database mass'
                        ' cut, greatly speed up the database loading and search but you can incur'
                        ' into empty boxes [default: None]')
    # parser.add_argument('--lf', type=int, help='OPTIMIZATION: Leafe size for'
                        # ' each node of the KDTree [default: %(default)d]',
                        # default=128)
    parser.add_argument('--force-light', action='store_true', help='DEBUG: '
                        'force the reading of the light catalogue first')


    args = parser.parse_args()

    if args.n < 1:
        raise ValueError('The value in -n must be a positive integer')

    if args.niter < 0:
        raise ValueError('The number of iterations in --niter must be a'
                         ' positive integer or 0 (for non overlapping search)')

    return args



def hebb_CLI():
    import numpy as np
    import os

    try:
        def_path = os.environ['HEBB_DB_PATH']
    except KeyError:
        raise ValueError('The environment variable HEBB_DB_PATH must be set.'
                         ' Check the "Database Setup" section at'
                         ' https://github.com/NegriAndrea/hebb')

    args = define_args()

    if args.l:
        print(listSnapshotTimes())
        return


    Mmax, fileNrMax, subNrMax = hebb(args.z_target, def_path,
                                     survey = (None if args.survey is
                                               None else tuple(args.survey)),
                                     NboxesOS = args.niter,
                                     L=args.L, M=args.M,
                                     force_light=args.force_light,
                                     NMassRank=args.n)
    if args.bb > 0:
        tvar = bootstrap_driver(Mmax, args.bb)
    else:
        tvar = None

    tSigmas, tRound, tErrors = buildPercentileTable(Mmax)
    prettyPrint(tSigmas, tRound, tErrors, write=args.t, tvar = tvar)


    if args.t:
        save_table(Mmax, fileNrMax, subNrMax, args)


    if args.plot:
        pretty_plot(Mmax)

def hebb_evolution_CLI():
    import numpy as np
    import os
    # mute the standard hebb logger
    logging.getLogger("hebb.logger_mine").setLevel(logging.INFO + 1)

    logger = logging.getLogger(__name__)

    try:
        def_path = os.environ['HEBB_DB_PATH']
    except KeyError:
        raise ValueError('The environment variable HEBB_DB_PATH must be set.'
                         ' Check the "Database Setup" section at'
                         ' https://github.com/NegriAndrea/hebb')

    description = (
            'Compute the N most massive dark matter haloes that you can find in a given '
            'survey with [z_min, z_max] and field-of-view with '
            'for all the snapshots available in Uchuu (2/h cGpc)^3 run.'
            )

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group( "Processing mode (required, mutually exclusive)")
    group_container = group.add_mutually_exclusive_group(required=True)
    group_container.add_argument('--survey', nargs=3, type=float,
                                 help='Survey min z, max z, FOV in arcmin^2'
                                 ' (used to compute box volume)',
                                 metavar=("z_min", "z_max", "fov"))
    group_container.add_argument('-L', type=float, help='Size of the box in cMpc,'
                                 ' alternative to the boxsize computation from the FOV'
                                 ' and z-depth of the survey')

    parser.add_argument('-n', type=int, help='Track the N most massive haloes'
                        ' in each box (1 tracks only the most massive) [default: %(default)d]', default=1)
    parser.add_argument('--force-light', action='store_true', help='DEBUG: '
                        'force the reading of the light catalogue first')


    args = parser.parse_args()

    if args.n < 1:
        raise ValueError('The value in -n must be a positive integer')


    snaps = listSnapshotTimes()

    tab_list = []
    tab_tSigma = []
    logger.info(f"EVOLUTION: starting computation for {len(snaps)} different z")

    for S, z in zip(snaps['SnapNr'].value, snaps['z'].value):

        # For each z try to use the largest mass cut possible: do a loop and if
        # we get an empty box divide the mass cut by a factor of 2
        M = 1e14
        for j in range(20):
            try:
                Mmax, fileNrMax, subNrMax = hebb(z, def_path,
                                                 survey = (None if args.survey is
                                                           None else tuple(args.survey)),
                                                 L=args.L, M=M,
                                                 force_light=args.force_light,
                                                 NMassRank=args.n)

            except EmptyBoxError:
                M/=2
                continue

            t = save_table(Mmax, fileNrMax, subNrMax, args)
            t['SnapNr'] = np.full(len(t), S, dtype=np.uint8)
            t['z'] = np.full(len(t), z, dtype=np.float32)
            tab_list.append(t)

            tSigmas, tRound, _ = buildPercentileTable(Mmax, units=True)
            tSigmas['SnapNr'] = np.full(len(tSigmas), S, dtype=np.uint8)
            tSigmas['z'] = np.full(len(tSigmas), z, dtype=np.float32)
            tab_tSigma.append(tSigmas)
            break

        logger.info(f"EVOLUTION: done for snapNr={S:2d} z={float(z):.2f}")

    t = vstack(tab_list)
    t.write('hebb_evolution.txt', format='ascii.ecsv', overwrite=True)

    tSigmas = vstack(tab_tSigma)
    for cname in tSigmas.colnames:
        tSigmas.rename_column(cname, cname.replace('σ', 'sig'))
    tSigmas.write('hebb_evolution_sigma_percentiles.txt', format='ascii.ecsv', overwrite=True)



def hebb_trace_CLI():

    description = ("Trace back galaxies found with hebb.")
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('hebbTable', type=str, help='Data table produced by hebb')
    parser.add_argument('z_target', type=float, help='Target redshift')
    parser.add_argument('treePath', type=str, help='Path of the merger tree file')


    parser.add_argument('-v', action='count', default=0,
                        help='verbosity level [%(default)d]')
    parser.add_argument('-s', action='store_true',
                        help='Use the code in serial, without the need of MPI')
    parser.add_argument('-f', action='store_true',
                        help='For each merger tree file, load the whole forest'
                        ' (best for when you have many galaxies in the same'
                        ' file)')

    args = parser.parse_args()

    finalT, mpirank = hebb_trace(args.hebbTable, args.treePath, args.z_target,
            args.v, args.s, readEntireForest = args.f)


    if mpirank == 0:
        finalT.write('hebb_traceback.txt', format='ascii.ecsv', overwrite=True)

def hebb_estimate_CLI():
    import numpy as np
    import os

    try:
        def_path = os.environ['HEBB_DB_PATH']
    except KeyError:
        raise ValueError('The environment variable HEBB_DB_PATH must be set.'
                         ' Check the "Database Setup" section at'
                         ' https://github.com/NegriAndrea/hebb')

    description = (
            'Estimate the smallest L at a certain z with a certain mass cut'
            )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('z_target', type=float, help='Number of sampling')
    parser.add_argument('L', type=float, help='Size of the box in cMpc')

    parser.add_argument('-v', action='count', default=0,
                        help='verbosity level [%(default)d]')

    parser.add_argument('-M', type=float, help='OPTIMIZATION: Database mass'
                        ' cut, greatly speed up the search but you can incour'
                        ' into empty boxes [default: None]')


    args = parser.parse_args()

    hebb_estimate(args.z_target, def_path, args.L, M=args.M, v=args.v)
