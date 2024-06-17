#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 20:53:54 2023

@author: nolan
"""

import numpy as np
import pandas as pd
import sys
import argparse

def main():
    """Convert PDS radar spectra to bare CSV

    Read a series of PDS-format radar spectra (per-run csv) and generate 
    simple CSVs containing each spectrum and the sdev-weighted sum.

    Returns
    -------
    None.

    """

    parser = argparse.ArgumentParser(
        description="convert PDS radar spectra to simple CSVs"
                    " Generates three output files:\n"
                    "spec.oc and spec.sc are simple csv tables with calibrated spectra\n"
                    "spec.sum contans the weighed sum")
    parser.add_argument('file', nargs='+',
                        help="Files to convert")

    parser.add_argument("-b", "--basename", action="store",
                        help="File stem to construct output files"
                        "Default is \"spec\" ",
                        default="spec")
    parser.add_argument("-n", "--noweight", action="store_true",
                        help="Don't weight sum by variance")
    parser.add_argument("-u", "--uncalibrated", action="store_true",
                        help="Don't apply calibration factors. Gives unweighted sum")
    args = parser.parse_args()





    wtsum1=0
    wtsum2=0
    first=True
    for i in args.file:
        try:
            d=pd.read_csv(i,header=None)
        except:
            print("Unable to parse csv file %s, quitting" %i, file=sys.stderr)    
            sys.exit(-1)
        try:
            sdevline=d[d[0]=='sdev'].index.tolist()[0]
        except IndexError:
            print("Couldn't find sdev in file %s. Quitting" % i, file=sys.stderr)
            sys.exit(-1)
        sdrow=d[:][sdevline:sdevline+1]
        sdev1=float(sdrow[1])
        sdev2=float(sdrow[2])

        if args.uncalibrated:
            fac1 = 1
            fac2 = 1
        else:
            if sdev1 <= 0 or sdev2 <= 0:
                print("Zero or negative sdev in file"
                      " %s, can't calibrate. Quitting." % i)
                sys.exit(-1)
            fac1 = sdev1
            fac2 = sdev2
        if args.noweight:
            wt1 = 1
            wt2 = 1
        else:
            wt1 = 1/fac1/fac1
            wt2 = 1/fac2/fac2
        wtsum1 += wt1
        wtsum2 += wt2
        try:
            datastart=d[d[0]=='# Data'].index.tolist()[0]+1
        except IndexError:
            print("Couldn't find data start marker in file"
                           " %s. Quitting" % i)
            sys.exit(-1)
        data=d[:][datastart:]
        pol1=data[1].to_numpy(dtype=float)
        pol2=data[2].to_numpy(dtype=float)
        freq=data[0].to_numpy(dtype=float)
        if first:
            sum1 = pol1*wt1*fac1
            sum2 = pol2*wt2*fac2
            arr1 = pol1*fac1
            arr2 = pol2*fac2
            first = False
            nf = freq.size
        else:
            sum1 += pol1*wt1*fac1
            sum2 += pol2*wt2*fac2
            arr1 = np.vstack((arr1,pol1*fac1))
            arr2 = np.vstack((arr2,pol2*fac2))
            if nf != freq.size:
                print(("Number of channels in file %s is %d but earlier "
                "files have %d") % (i, freq.size,nf))
                sys.exit(-1)
    sum1 = sum1 / wtsum1
    sum2 = sum2 / wtsum2
    df1=pd.DataFrame(data=np.vstack((freq,arr1)).transpose(),
                     columns=['Frequency']+args.file)
    df2=pd.DataFrame(data=np.vstack((freq,arr2)).transpose(),
                     columns=['Frequency']+args.file)
    df3=pd.DataFrame(data=np.vstack((freq,sum1,sum2)).transpose(),
                     columns=['Frequency','Pol1','Pol2'])
    df1.to_csv(args.basename+".oc.csv",index=False)
    df2.to_csv(args.basename+".sc.csv",index=False)
    df3.to_csv(args.basename+".sum.csv",index=False)

if __name__ == "__main__":
    main()
