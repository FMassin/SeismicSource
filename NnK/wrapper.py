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

# def ggg(...):
#         """
#         Run a ...
#
#         :param type: String that specifies which trigger is applied (e.g.
#             ``'recstalta'``).
#         :param options: Necessary keyword arguments for the respective
#             trigger.
#
#         .. note::
#
#             The raw data is not accessible anymore afterwards.
#
#         .. rubric:: _`Supported Trigger`
#
#         ``'classicstalta'``
#             Computes the classic STA/LTA characteristic function (uses
#             :func:`obspy.signal.trigger.classicSTALTA`).
#
#         .. rubric:: Example
#
#         >>> ss.ssss('sss', ss=1, sss=4)  # doctest: +ELLIPSIS
#         <...aaa...>
#         >>> aaa.aaa()  # aaa
#
#         .. plot::
#
#             from ggg import ggg
#             gg = ggg()
#             gg.ggggg("ggggg", ggggg=3456)
#             gg.ggg()
#         """

import os

def readallchannels(dataset, operation='eventdir'):
    """
    wrapps obspy.core.stream.read so various seismic file 
    formats can be read in one pass.

    Assumes data files are organized in event directories.
    This can be improved
    
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

    
    # Incremental slurps 
    md=-1
    for name in files: # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
        #print "Reading",name,"..."
        with open (name, "r") as myfile:
            dataset.extend(myfile.readlines()) # get ride of \n !!!
        
        if operation is 'relative':
            pathname=os.path.dirname(name)
            for d, data in enumerate(dataset):
                if d >= md:
                    dataset[d]=pathname+'/'+data.strip()
                    md=d

    ## Stream processing version
    ## the stream version seems less efficient than the slurp version
    ## probably due to the incremental increase of list dataset
    ## Also something is wrong in the returned dataset
    #
    # for fname in files:
    #     print "Reading",fname,"..."
    #     if operation is None :
    #         with open(fname, 'r+') as f:
    #             for line in f:
    #                 dataset.extend(line)

    #     if operation is 'relative':
    #         pathname=os.path.dirname(fname)
    #         with open(fname, 'r+') as f:
    #             for line in f:                    
    #                 dataset.extend(pathname+'/'+line)    

    return dataset

def randomsample(dataset,Nsample,searchregex=None, test=None) :
    """
    Wrapps random.sample so list elements are randomly
    resampled to the given dimension (2nd arg).

    Optionally filter the sample with regular expression.

    Optionally filter the sample with only existing full
    file names in elements.

    
    :type dataset: list
    :param dataset: files (patterns) to resample.
    :type Nsample: Variable
    :param Nsample: total number in returned sample.
    :type searchregex: String
    :param searchregex: regular expression for filtering.
    :type test: String
    :param test: None or 'existingonly' to test file existence.
    :rtype: list
    :return: contents of files.
    """
    # 1) Test files and search patterns.
    # 2) Read files and add to returned list.
    # 3) Concatenate list and catalog paths if asked.
    
    import random
    import re

    # Usual instance of random.sample, if no optinal test are 
    # requested, the wrapper is not usefull at all
    tempdataset = random.sample(dataset, Nsample)

    if (test is 'existingonly') or (searchregex):

        for e, eventfile in enumerate(tempdataset):  
            
            test1 = True
            if test is 'existingonly':
                test1 = os.path.exists(eventfile)

            test2 = re.search(searchregex, eventfile)
            
            while (not test1) or (not test2):
                
                eventfile = random.sample(dataset, 1)
                eventfile=eventfile[0]
                
                test1 = True
                if test is 'existingonly':
                    test1 = os.path.exists(eventfile)

                test2 = re.search(searchregex, eventfile)
            
            tempdataset[e]=eventfile

    return tempdataset

