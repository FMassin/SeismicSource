# -*- coding: utf-8 -*-
"""
wrapper - Module for wrapping of existing functions 
======================================================
This module ...

.. figure:: /_images/Event.png

.. note::

    For ...

:copyright:
    The ...
:license:
    ...
"""

import os

def readallchannels(dataset):
    """
    wrapps obspy.core.stream.read so various seismic file 
    formats can be read in one pass.

    :type dataset: list
    :param dataset: files to read.
    :rtype: class obspy.core.stream
    :return: class of data streams.
    """
    # 1) Iterate over full files names
    # 2) Get all available channel of each file 
    # 3) Read individually in stream.

    from obspy.core.stream import Stream, read
    import glob
    import sys

    # Reading the waveforms
    m=0
    eventset=[]
    waveformset = Stream()
    for e, eventfile in enumerate(dataset):  

        eventfile = eventfile.strip()
        waveformset += read(eventfile)
        eventset.append(e)
        (root,ext)=os.path.splitext(eventfile)
        channelset = glob.glob(os.path.dirname(eventfile)+'/*'+waveformset[-1].stats.station+'*'+ext)
        
        for c, channelfile in enumerate(channelset):    
            if not eventfile == channelfile :
                eventset.append(e)
                waveformset += read(channelfile) # file with uncorrect metadata may be ignored (ex: YHHN in station metatdata)
                
                #df = waveformset[-1].stats.sampling_rate
                #cft = classicSTALTA(waveformset[-1].data, int(5 * df), int(10 * df))
                #plotTrigger(waveformset[-1], cft, 1.5, 0.5)

    return waveformset

def readfullfilenames(paths, operation=None):
    """
    Wrapps readlines so several text files could be 
    loaded in one pass.

    Reads full file names and paths in given catalog 
    files.

    Optionally stitch catalog content with the path of
    the given catalog.

    
    :type dataset: list
    :param dataset: files (patterns) to read.
    :rtype: list
    :return: contents of files.
    """
    # 1) Test files and search patterns.
    # 2) Read files and add to returned list.
    # 3) Concatenate list and catalog paths if asked.
    
    import glob


    files=[]
    dataset=[]
    for p in paths:
        files.extend(glob.glob(p))

    md=-1
    for name in files: # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
        #print "Reading",name,"..."
        with open (name, "r") as myfile:
            dataset.extend(myfile.readlines()) # get ride of \n !!!
        
        if operation is 'relative':
            pathname=os.path.dirname(name)
            for d, data in enumerate(dataset):
                if d >= md:
                    dataset[d]=pathname+'/'+data
                    md=d

    return dataset

def randomsample(dataset,Nsample,searchregex=".*", test=None) :
    """
    Wrapps random.sample so list elements are randomly
    resampled to the given dimension (2nd arg).

    Optionally filter the sample with regular expression.

    Optionnally filter the sample with only existing full
    file names in elements.

    
    :type dataset: list
    :param dataset: files (patterns) to resample.
    :type Nsample: Variable
    :param Nsample: total number in returned sample.
    :type searchregex: String
    :param searchregex: regular expression for filtering.

    :rtype: list
    :return: contents of files.
    """
    # 1) Test files and search patterns.
    # 2) Read files and add to returned list.
    # 3) Concatenate list and catalog paths if asked.
    
    import random
    import re

    tempdataset = random.sample(dataset, Nsample)
    for e, eventfile in enumerate(tempdataset):  
        
        test1 = True
        if test is 'existingonly':
            test1 = os.path.exists(eventfile.strip())

        test2 = re.search(searchregex, eventfile.strip())
        
        while (not test1) or (not test2):
            
            eventfile = random.sample(dataset, 1)
            eventfile=eventfile[0]
            
            test1 = True
            if test is 'existingonly':
                test1 = os.path.exists(eventfile.strip())

            test2 = re.search(searchregex, eventfile.strip())
        
        tempdataset[e]=eventfile
        #print eventfile.strip()

    return tempdataset

