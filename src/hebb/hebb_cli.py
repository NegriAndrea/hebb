#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# import numpy as np
# import matplotlib.pyplot as plt
# import h5py
# from pathlib import PurePath
from .hebb_core import hebb

def hebb_CLI():
    import numpy as np
    import argparse

    description = (
            'Compute the heaviest dark matter halo that you can find in a given '
            'survey with [z_min, z_max] and field-of-view by performing a '
            'non-parametric block bootstrap over the Uchuu (2/h cGpc)^3 run.'
            )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('Nboxes', type=int, help='Number of boxes for bootstrap')
    parser.add_argument('z_target', type=float, help='Number of sampling')
    parser.add_argument('z_min', type=float, help='')
    parser.add_argument('z_max', type=float, help='')
    parser.add_argument('fov', type=float, help='FOV in arcmin^2')

    parser.add_argument('-v', action='count', default=0,
                        help='verbosity level [%(default)d]')
    parser.add_argument('-p', type=str, help='Path of the catalogue file'
                        ' [default: %(default)s]', default='~/hebb_data')
    parser.add_argument('-t', action='store_true', help='Create a table with'
                        ' the sampled haloes')
    parser.add_argument('--plot', action='store_true', help='show a plot of the M200 distribution')

    parser.add_argument('-L', type=float, help='Size of the box in cMpc,'
                        ' overrides the boxsize computation from the FOV')
    parser.add_argument('-M', type=float, help='OPTIMIZATION: Database mass'
                        ' cut, greatly speed up the search but you can incour'
                        ' into empty boxes [default: None]')
    parser.add_argument('--lf', type=int, help='OPTIMIZATION: Leafe size for'
                        ' each node of the KDTree [default: %(default)d]',
                        default=128)


    args = parser.parse_args()

    Mmax, fileNrMax, subNrMax = hebb(args.z_target, args.Nboxes,
                                     args.p, args.z_min, args.z_max, args.fov,
                                     L=args.L, M=args.M, v=args.v, leafsize =
                                     args.lf)

    M16, M50, M86 = np.quantile(10.**Mmax, [0.16, 0.5, 0.86])
    print(f"log10(M200) = {np.log10(M50)} _-{np.log(10)*(M50-M16)/M50} ^+{np.log(10)*(M86-M50)/M50}")
    print(f"log10(M200) = {np.log10(M50):.2f} _-{np.log(10)*(M50-M16)/M50:.2f} ^+{np.log(10)*(M86-M50)/M50:.2f}")

    if args.plot:
        import matplotlib.pyplot as plt
        plt.hist(Mmax, 50)
        plt.show()

    if args.table:
        from astropy.table import Table
        t=Table([Mmax, fileNrMax, subNrMax], names=['M200', 'fileNr', 'subNr'])
        t.write('table_max.txt', format='ascii.ecsv', overwrite=True)

