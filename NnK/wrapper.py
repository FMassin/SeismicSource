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


def read(a):
    """
    Calculate the dimensions of all data in stream.

    Given a stream (obspy.core.stream) calculate the minimum dimensions 
    of the array that contents all data from all traces.

    This method is written in pure Python and gets slow as soon as there
    are more then ... in ...  --- normally
    this does not happen.

    :type a: ObsPy :class:`~obspy.core.stream`
    :param a: datastream of e.g. seismogrammes.
    :rtype: array
    :return: array of int corresponding to the dimensions of stream.
    """
    # 1) Does this
    # 2) Then does 
    #    that