#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Adjust a series of sband images to have the same eph_row.

Default is to choose the first image.
It will error if asked to shift by more than one whole images, but does not
otherwise know that the shift is reasonable.

Created on Fri Mar 13 16:55:02 2020
@author: nolan
"""

import argparse
from astropy.io import fits
import os
import shutil
import scipy.ndimage
import numpy


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

    if args.EPHROW:
        erow = args.EPHROW
    else:
        ff = args.filename[0]
        hdl = fits.open(ff, mode='readonly')
        if args.goal and 'GOAL_ROW' in hdl[0].header:
            erow = hdl[0].header['GOAL_ROW']
            print(f"Aligning to GOAL_ROW of {erow}.")
        elif 'EPH_ROW' in hdl[0].header:
            erow = hdl[0].header['EPH_ROW']
            print(f"Aligning to EPH_ROW of {erow}.")
        else:
            raise KeyError('No EPH_ROW found or specified')

    for f in args.filename:
        if args.inplace:
            outfile = ''
        else:
            base, ext = os.path.splitext(f)
            outfile = base + args.suffix + ext

        sbalign(f, outfile, erow, args.inplace, args.boundary, args.cval)


def sbalign(infile, outfile, erow, inplace=False, boundary='wrap', cval=0):
    """Sbalign one file.

    Inputs
    ------
    input file name
    output file name (ignoded if in-place)
    inplace Boolean
    erow - Epehemris row to align to
    boundary condition
    constant value for boundary = 'const'


    Returns: None

    Has side effect of adding descriptive keywords even if no shift done
    """
    if inplace:
        hdl = fits.open(infile, 'update')
    else:
        shutil.copyfile(infile, outfile)
        hdl = fits.open(outfile, 'update')

        try:
            ferow = hdl[0].header['EPH_ROW']
        except KeyError:
            raise KeyError(f'No EPH_ROW found in file {hdl.info()}')
        ysize = hdl[0].header['NAXIS2']
        shift = erow - ferow
        shift = round(shift)
        newerow = ferow + shift  # Leave fraction in eph_row
        if abs(shift) > ysize-1:
            hdl.close()
            raise UserWarning(f'Would have to shift {infile} by {shift}'
                              ' pixels but it is only {ysize} pixels tall.\n'
                              'EPH_ROW is {ferow}, shifting to {newerow},'
                              ' skipping')
        hdl[0].header['EPH_ROW'] = newerow
        hdl[0].header['CRPIX2'] = newerow + 1
        hdl[0].header['GOAL_ROW'] = (erow, 'Goal row in sbalign')

        if not ('ORIG_ROW' in hdl[0].header):  # There is one original row
            hdl[0].header['ORIG_ROW'] = ferow
        hdl[0].header['HISTORY'] = (f'Shift of {shift} applied in sbalign')
        data = hdl[0].data
        print(f'Shifting {infile} by {shift} pixels.')
        if 'wrap' == boundary:  # scipy.ndimage has a known bug with wrap
            sdata = numpy.roll(data, shift, axis=0)
        elif 'const' == boundary:
            sdata = scipy.ndimage.shift(data, (shift, 0), order=0,
                                        mode='constant', cval=cval)
        else:
            sdata = scipy.ndimage.shift(data, (shift, 0), order=0,
                                        mode=boundary)
        hdl[0].data = sdata
        hdl.close()


if __name__ == '__main__':
    main()
