#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from .hebb_core import hebb, hebb_estimate
from .hebb_trace import hebb_trace
from .utils import pretty_print
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[
        # logging.FileHandler("hebb_run.log"),
        logging.StreamHandler()
    ]
)



def hebb_CLI():
    import numpy as np
    import argparse
    import os

    try:
        def_path = os.environ['HEBB_DB_PATH']
    except KeyError:
        raise ValueError('The environment variable HEBB_DB_PATH must be set.'
                         ' Check the "Database Setup" section at'
                         ' https://github.com/NegriAndrea/hebb')

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

    if args.survey is None:
        z1 = None
        z2 = None
        fov = None
    else:
        z1  = args.survey[0]
        z2  = args.survey[1]
        fov = args.survey[2]

    if args.n < 1:
        raise ValueError('The value in -n needs to be a positive integer')

    Mmax, fileNrMax, subNrMax = hebb(args.z_target, args.Nboxes,
                                     def_path, z1, z2, fov,
                                     L=args.L, M=args.M, leafsize =
                                     args.lf, force_light=args.force_light,
                                     nn=args.n)
    pretty_print(Mmax)


    if args.t:
        import astropy.units as u
        from astropy.table import Table
        # unroll them
        ordering = np.repeat(np.arange(Mmax.shape[0],
                                       dtype=np.min_scalar_type(args.n)),
                             Mmax.shape[1])

        # alter shape instead to use np.reshape to force an error
        # in case the code attempts to copy the arrays (should not happen)
        fileNrMax.shape = fileNrMax.size
        subNrMax.shape = subNrMax.size
        ordering.shape = ordering.size
        # take a view to have it 1d and not to mess with the plot below
        # no copy is involved
        Mmax1d = Mmax.view()
        Mmax1d.shape = Mmax.size

        t=Table([Mmax1d, ordering, fileNrMax, subNrMax],
                names=['M200', 'MassRank', 'fileNr', 'subNr'])

        t['M200'].unit = u.dex(u.Msun)
        t['MassRank'].description = ('Rank of the halo in mass sorting'
                            f" ([0..{args.n-1}], 0 is the most massive)")

        t.meta = {'Description':'Hebb result table','Nboxes':args.Nboxes,
                  'z_target':args.z_target, '-n':args.n,
                  '--survey':args.survey, '-L':args.L, '-M': args.M}

        t.sort(['MassRank', 'fileNr', 'subNr'])
        t.write('table_max.txt', format='ascii.ecsv', overwrite=True)


    if args.plot:
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
