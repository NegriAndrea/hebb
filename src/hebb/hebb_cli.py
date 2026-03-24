#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .hebb_core import hebb, hebb_estimate
from .hebb_trace import hebb_trace
from .utils import pretty_print, save_table, pretty_plot
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[
        # logging.FileHandler("hebb_run.log"),
        logging.StreamHandler()
    ]
)

def define_args() -> argparse.Namespace:
    """
    Define the command line interface and options of `hebb`.

    """

    description = (
            'Compute the N most massive dark matter haloes that you can find in a given '
            'survey with [z_min, z_max] and field-of-view by performing a '
            'non-parametric block bootstrap over the Uchuu (2/h cGpc)^3 run.'
            )

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('Nboxes', type=int, help='Number of boxes for bootstrap')
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
    parser.add_argument('-t', action='store_true', help='Create a table with'
                        ' the sampled haloes')
    parser.add_argument('--plot', action='store_true', help='Show a plot of the M200 distribution')

    parser.add_argument('-M', type=float, help='OPTIMIZATION: Database mass'
                        ' cut, greatly speed up the database loading and search but you can incur'
                        ' into empty boxes [default: None]')
    parser.add_argument('--lf', type=int, help='OPTIMIZATION: Leafe size for'
                        ' each node of the KDTree [default: %(default)d]',
                        default=128)
    parser.add_argument('--force-light', action='store_true', help='DEBUG: '
                        'force the reading of the light catalogue first')


    args = parser.parse_args()

    if args.n < 1:
        raise ValueError('The value in -n needs to be a positive integer')

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


    Mmax, fileNrMax, subNrMax = hebb(args.z_target, args.Nboxes, def_path,
                                     survey = (None if args.survey is
                                               None else tuple(args.survey)),
                                     L=args.L, M=args.M, leafsize =
                                     args.lf, force_light=args.force_light,
                                     NMassRank=args.n)
    pretty_print(Mmax, write=args.t)


    if args.t:
        save_table(Mmax, fileNrMax, subNrMax, args)


    if args.plot:
        pretty_plot(Mmax)

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

    args = parser.parse_args()

    finalT, mpirank = hebb_trace(args.hebbTable, args.treePath, args.z_target,
            args.v, args.s)


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
