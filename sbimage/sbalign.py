#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sbalign will adjust a series of sband images to have the same eph_row.

Default is to choose the first image.

Created on Fri Mar 13 16:55:02 2020
@author: nolan
"""

import argparse
from astropy.io import fits
import sys
import os
import shutil
import scipy.ndimage


def main():
    """Run sbalign.

    Returns
    -------
    None.

    """
    parser = argparse.ArgumentParser(
        description="Align images on EPH_ROW")
    parser.add_argument('filename', nargs='+',
                        help="Files to align")

    parser.add_argument("-E", "--EPHROW", action="store", type=float,
                        help="EPH_ROW to apply to all images.",
                        )
    parser.add_argument('-g', '--goal', action='store_true',
                        help='Use GOAL_ROW in first image rather than EPH_ROW'
                        ' to match a previous use of sbalign')
    parser.add_argument("-b", "--boundary", choices=['wrap', 'nearest',
                                                     'const'],
                        default='wrap',
                        help='how to treat edges')
    parser.add_argument("-c", "--cval", action="store", default=0, type=float,
                        help='Constant to use for constant padding, default 0')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--inplace", action="store_true",
                       help='edit files in-places rather than appending .s')
    group.add_argument('-s', '--suffix', action='store', default='.s',
                       help='suffix to add to new file name')
    args = parser.parse_args()

    first = True
    for f in args.filename:
        if args.inplace:
            hdl = fits.open(f, mode='update')
        else:
            base, ext = os.path.splitext(f)
            newfile = base + args.suffix + ext
            shutil.copyfile(f, newfile)
            hdl = fits.open(newfile, 'update')

        if first:
            first = False
            if args.EPHROW:
                erow = args.EPHROW
            else:
                nogoal = True
                if args.goal:
                    nogoal = False
                    try:
                        erow = hdl[0].header['GOAL_ROW']
                    except KeyError:
                        nogoal = True
                if nogoal:
                    try:
                        erow = hdl[0].header['EPH_ROW']
                    except KeyError:
                        exit('No EPH_ROW found or specified')
            print(f'Aligning to new EPH_ROW of {erow}')

        try:
            ferow = hdl[0].header['EPH_ROW']
        except KeyError:
            sys.exit(f'No EPH_ROW found in file {hdl.info()}')
        ysize = hdl[0].header['NAXIS2']
        shift = erow - ferow
        shift = round(shift)
        newerow = ferow + shift  # Leave fraction in eph_row
        if abs(shift) > ysize-1:
            print(f'Would have to shift {f.name} by {shift} pixels but it is',
                  ' only {ysize} pixels tall.')
            print(f'EPH_ROW is {ferow}, shifting to {newerow}, skipping')
            hdl.close()
            next

        hdl[0].header['EPH_ROW'] = newerow
        hdl[0].header['CRPIX2'] = newerow + 1
        hdl[0].header['GOAL_ROW'] = (erow, 'Goal row in sbalign')
        hdl[0].header['HISTORY'] = (f'Shift of {shift} applied in sbalign')
        data = hdl[0].data
        print(f'Shifting {f} by {shift} pixels.')
        if args.boundary == 'const':
            sdata = scipy.ndimage.shift(data, (shift, 0), mode='constant',
                                        cval=args.cval)
        else:
            sdata = scipy.ndimage.shift(data, (shift, 0), mode=args.boundary)
        hdl[0].data = sdata
        hdl.close()


if __name__ == '__main__':
    main()
