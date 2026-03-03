#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .hebb_core import hebb, hebb_estimate
from .hebb_trace import hebb_trace

def hebb_CLI():
    import numpy as np
    import argparse
    import os

    try:
        def_path = os.environ['HEBB_DB_PATH']
    except KeyError:
        raise ValueError('The environment variable HEBB_DB_PATH must be set')

    description = (
            'Compute the heaviest dark matter halo that you can find in a given '
            'survey with [z_min, z_max] and field-of-view by performing a '
            'non-parametric block bootstrap over the Uchuu (2/h cGpc)^3 run.'
            )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('Nboxes', type=int, help='Number of boxes for bootstrap')
    parser.add_argument('z_target', type=float, help='Number of sampling')
    # parser.add_argument('z_min', type=float, help='Survey min z (used to compute box volume)')
    # parser.add_argument('z_max', type=float, help='Survey max z (used to compute box volume)')
    # parser.add_argument('fov', type=float, help='FOV in arcmin^2 (used to compute box volume)')

    # requiredNamed = parser.add_argument_group('required named arguments')
    group = parser.add_argument_group( "Processing mode (required, mutually exclusive)")
    group_container = group.add_mutually_exclusive_group(required=True)
    group_container.add_argument('--survey', nargs=3, type=float,
                                 help='Survey min z, max z, fov in arcmin^2'
                                 ' (used to compute box volume)',
                                 metavar=("z_min", "z_max", "fov"))
    group_container.add_argument('-L', type=float, help='Size of the box in cMpc,'
                        ' overrides the boxsize computation from the FOV')

    parser.add_argument('-v', action='count', default=0,
                        help='verbosity level [%(default)d]')
    # parser.add_argument('-p', type=str, help='Path of the catalogue file'
                        # ' [default: %(default)s]', default=def_path)
    parser.add_argument('-t', action='store_true', help='Create a table with'
                        ' the sampled haloes')
    parser.add_argument('--plot', action='store_true', help='show a plot of the M200 distribution')

    parser.add_argument('-M', type=float, help='OPTIMIZATION: Database mass'
                        ' cut, greatly speed up the search but you can incour'
                        ' into empty boxes [default: None]')
    parser.add_argument('--lf', type=int, help='OPTIMIZATION: Leafe size for'
                        ' each node of the KDTree [default: %(default)d]',
                        default=128)


    args = parser.parse_args()

    if args.survey is None:
        z1 = None
        z2 = None
        fov = None
    else:
        z1  = args.survey[0]
        z2  = args.survey[1]
        fov = args.survey[2]

    Mmax, fileNrMax, subNrMax = hebb(args.z_target, args.Nboxes,
                                     def_path, z1, z2, fov,
                                     L=args.L, M=args.M, v=args.v, leafsize =
                                     args.lf)

    M16, M50, M86 = np.quantile(10.**Mmax.astype(np.float64), [0.16, 0.5, 0.86])
    print(f"log10(M200) = {np.log10(M50)} _-{(M50-M16)/M50/np.log(10)} ^+{(M86-M50)/M50/np.log(10)}")
    print(f"log10(M200) = {np.log10(M50):.2f} _-{(M50-M16)/M50/np.log(10):.2f} ^+{(M86-M50)/M50/np.log(10):.2f}")

    if args.plot:
        import matplotlib.pyplot as plt
        plt.hist(Mmax, 50)
        plt.show()

    if args.t:
        from astropy.table import Table
        t=Table([Mmax, fileNrMax, subNrMax], names=['M200', 'fileNr', 'subNr'])
        t.sort(['fileNr', 'subNr'])
        t.write('table_max.txt', format='ascii.ecsv', overwrite=True)

def hebb_trace_CLI():
    import argparse

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
    import argparse
    import os

    try:
        def_path = os.environ['HEBB_DB_PATH']
    except KeyError:
        raise ValueError('The environment variable HEBB_DB_PATH must be set')

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
