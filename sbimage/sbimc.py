#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Combine Sbandin images.

Combine sband images into sums, keeping track of sband keywords.
Created on Mon Mar 16 22:36:58 2020

@author: nolan
"""

import argparse
from astropy.io import fits
import os
import numpy


def sbimc(files, outfile):
    """Sum sband images.

    Parameters
    ----------
    files : list of strings
        DESCRIPTION.
    outfile : string
        DESCRIPTION.

    Returns
    -------
    None.

    """
    exposure = 0.0
    jdstart = 999999999.
    jdend = -999999999.
    looks = 0
    tsyssum = 0.
    sigsum = 0.
    jdmeansum = 0

    if os.file.exists(outfile):
        raise(FileExistsError(f'File Exists: {outfile}'))

    first = True
    nels = len(files)
    weights = numpy.ones(nels)  # in case we want to implement
    sw = sum(weights)
    for f, w in zip(files, weights):
        hdu = fits.open(f, 'readonly')
        h = hdu[0].header
        exposure += h['EPOSURE']
        jdmeansum += h['JDMEAN'] * w
        looks += h['LOOKS']
        jdstart = min(jdstart, h['JDSTART'])
        jdend = max(jdend, h['JDEND'])
        tsyssum += h['TSYSSIGS']
        sigsum += h['SIGCNTS']
        if first:
            d = f[0].data * w
            saveheader = h
        else:
            d = d + f[0].data * w
        hdu.close()

    hdu = fits.open(outfile, 'update')
    hdu[0].header = saveheader
    hdu[0].data = d / sw
    h = hdu[0].header
    h['EXPOSURE'] = exposure
    h['JDMEAN'] = jdmeansum / sw
    h['LOOKS'] = looks
    h['JDSTART'] = jdstart
    h['JDEND'] = jdend
    h['TSYSSIGS'] = tsyssum / sw
    h['SIGCNTS'] = sigsum / sw
    hdu.close()



def main():
    """


    Returns
    -------
    None.

    """
    parser = argparse.ArgumentParser(
        description="Combine Images")
    parser.add_argument('filename', nargs='+',
                        help="Files to combine")

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
