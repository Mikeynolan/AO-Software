#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 22:06:55 2023

@author: nolan
"""
import numpy as np
import pandas as pd

def readcsvspec(files,calibrate=True):
    """Ingest PDS Radar Spectra

    This routine reads in a series of PDS Radar spectra CSV files and parses
    the data and the metadata.
    
    Parameters
    ----------
    files : list of csv filenames
        List of one or more csv files in the PDS radar spectra format. 
    
    
    calibrate : bool (default True)
        By default the spectra are returned as calibrated cross-section by
        multiplying each spectrum by the sdev tag.
        if calbrate is False they will remain normalized to unit standard 
        deviation.
    

    Returns
    -------
    numpy array of frequency axes. dimensions are [numspec, nfreq]
    numpy array of spectra: dimensions are [numspec,2 pol, nfreq]
    array of dicts containing the PDS keywords in each csv file
    array of dicts containing the "tags" for each spectrum (two values each)
    array of dicts containing the "extra tags" for each spectrum

    """

    first=True
    for i in files:  
        d=pd.read_csv(i,header=None)
        #Can't search through NaNs, but removing causes other trouble
        ds=d.fillna('')
        # Gather up the markers, fail if not all found
        try:
            keysstart=d[ds[0].str.contains(r'(?i)# *Keywords')].index.tolist()[0]+1
        except IndexError:
            print("Couldn't find keywordstart marker in file %s." % i)
        try:    
            tagsstart=d[ds[0].str.contains(r'(?i)Tags')].index.tolist()[0]+1
        except IndexError:
            print("Couldn't find tags start marker in file %s." % i)
        try:
            extrasstart=d[ds[0].str.contains(r'(?i)Extra *Tags')].index.tolist()[0]+1
        except IndexError:
            print("Couldn't find Extra Tags start marker in file %s." % i)
        try:
            coldefstart=d[ds[0].str.contains(r'(?i)# *Column *Definitions')].index.tolist()[0]+1
        except IndexError:
            print("Couldn't find coldefs start marker in file %s." % i)
        try:    
            datastart=d[ds[0].str.contains(r'(?i)# *Data')].index.tolist()[0]+1
        except IndexError:
            print("Couldn't find data start marker in file %s." % i)
# Most ranges go from start to end-2, avoiding marker
        keywords=d[:][keysstart:tagsstart-1]
        tags=d[:][tagsstart:extrasstart-1]
        extras=d[:][extrasstart:coldefstart-1]
        data=d[:][datastart:].to_numpy(dtype=float)
        freq=data[:,0]
        spec = data[:,1:3].transpose()
        
        kd = dict(keywords.values[:,0:2].tolist())
        td = dict(zip(tags.values[:,0].tolist(),tags.values[:,1:3].tolist()))
        ed = dict(zip(extras.values[:,0].tolist(),extras.values[:,1].tolist()))
        if calibrate:
            if 'sdev' in td:
                spec[0,:] *= float(td['sdev'][0])
                spec[1,:] *= float(td['sdev'][1])
            else:
                print("Calibration requested but no sdev in file %s, quitting"%i)
                exit(-1)
        if first:
            freqa=np.array(freq,ndmin=2)
            speca=np.array(spec,ndmin=3)
            keysa=[kd]
            tagsa=[td]
            extsa=[ed]
            first=False
        else:
            freqa=np.append(freqa,[freq],0)
            speca=np.append(speca,[spec],0)
            keysa.append(kd)
            tagsa.append(td)
            extsa.append(ed)
    return(freqa,speca,keysa,tagsa,extsa)
