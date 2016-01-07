# -*- coding: utf-8 -*-
"""
trigger - Module to improve obspy.signal.trigger
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

def correlatecomponents(a, scales=None, operation='cecm'):
    """
    correlatecomponents performs correlation between the 
    components of the same seismometers by rolling 
    operation.

    :type a: ObsPy :class:`~obspy.core.stream`
    :param a: datastream of e.g. seismogrammes.
    :type scales: vector
    :param scales: scale(s) of correlation operation.
    :type operation: string
    :param operation: type of operation.
    :rtype: array
    :return: array of root mean square time series, scale 
    	in r.stats._format (samples scale unit).

    .. note::

        The raw data is not accessible anymore afterwards.

    .. rubric:: _`Supported Trigger`

    ``'classicstalta'``
        Computes the classic STA/LTA characteristic function (uses
        :func:`obspy.signal.trigger.classicSTALTA`).

    .. rubric:: Example

    >>> ss.ssss('sss', ss=1, sss=4)  # doctest: +ELLIPSIS
    <...aaa...>
    >>> aaa.aaa()  # aaa

    .. plot::

        from ggg import ggg
        gg = ggg()
        gg.ggggg("ggggg", ggggg=3456)
        gg.ggg()
    """
    # 1) Calculate moving rms and average, see moving()
    # 2) ...
    # 3) ... 

    import copy
    from obspy.core.stream import Stream
    import numpy as np
    import re

    # Initialize multiscale if undefined
    (tmax,nmax) = streamdatadim(a)
    if scales is None:
        scales = [2**i for i in range(3,999) if 2**i < (nmax - 2**i)]
        scales = np.require(scales, dtype=np.int) 
    
    # Initialize results at the minimal size
    correlations = np.zeros(( tmax, len(scales), nmax )) 