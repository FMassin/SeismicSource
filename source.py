# -*- coding: utf-8 -*-
"""
source - Module to model seismic sources
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

def rpgen(strike, dip, rake, gamma, sigma, TKO, AZM) : 
    """
    rpgen(strike,dip,rake,gamma,sigma, TKO, AZM) calculates P-wave, S-wave,
    SH-wave and SV-wave radiation pattern unp.sing shear-tensile source model
    presented in [see references 1, 2, 3 for details]. All input angles 
    (strike, dip, rake of the fault, tensile angle gamma, takeoff angle 
    TKO and azimuth from the source to the observation point AZM) should 
    be in degrees. The takeoff angle is measure from bottom. The function 
    returns matrices of the same size as input TKO and AZM matrices. 

    :type strike, dip, rake: ObsPy : variables
    :param strike, dip, rake: fault plane parameters (degrees).
    :type gamma: variable
    :param gamma: tensile angle in degrees (0 degrees for pure shear faulting, 
        90 degrees for pure tensile opening).
    :type sigma: variable
    :param sigma:  Poisson's ratio.
    :type TKO: variable
    :param TKO: matrix of takeoff angles for which to calculate the 
        corresponding radiation pattern coefficients (degrees, the takeoff 
        angles are measured from bottom).
    :type AZM: variable
    :param AZM: matrix of corresponding azimuths (in degrees) for which the 
        radiation pattern coefficients should be calculated.
    :rtype Gp, Gs, Gsh, Gsv: vectors
    :return Gp, Gs, Gsh, Gsv: P-wave, S-wave, SH-wave, and SV-wave radiation 
         pattern coefficients calculated for corresponding takeoff angles 
         and azimuths specified in TKO and AZM matrices.

    .. seealso::

        [1] Kwiatek, G. (2013/09/15). Radiation pattern from shear-tensile 
         seismic source. Revision: 1.3. http://www.npworks.com/matlabcentral/fileexchange/43524-radiation-pattern-from-shear-tensile-seismic-source

        [2] Kwiatek, G. and Y. Ben-Zion (2013). Assessment of P and S wave 
         energy radiated from very small shear-tensile seismic events in 
         a deep South African mine. J. Geophys. Res. 118, 3630-3641, 
         doi: 10.1002/jgrb.50274

        [3] Ou, G.-B., 2008, Seismological Studies for Tensile Faults. 
         Terrestrial, Atmospheric and Oceanic Sciences 19, 463.

        [4] Vavryèuk, V., 2001. Inversion for parameters of tensile 
         earthquakes.” J. Geophys. Res. 106 (B8): 16339–16355. 
         doi: 10.1029/2001JB000372.

    .. rubric:: Example

        import numpy as np
        
        AZIMUTH = np.arange(0,360,1)*np.pi/180
        TAKEOFF = np.arange(0,180,1)*np.pi/180
        [AZIMUTH,TAKEOFF] = np.meshgrid(AZIMUTH,TAKEOFF)

        strike = 50
        dip = 60
        rake = -90
        sigma = 0.25
        gamma = 0
        [GP, GS, GSH, GSV] = source.rpgen(strike,dip,rake,gamma,sigma, TAKEOFF*180/np.pi, AZIMUTH*180/np.pi)
               

    .. plot::

        import matplotlib
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        scale = np.abs(GP)
        YP = scale*np.cos(AZIMUTH) * np.sin(TAKEOFF)
        XP = scale*np.sin(AZIMUTH) * np.sin(TAKEOFF)
        ZP = scale*-np.cos(TAKEOFF)

        # initializing the colormap machinery
        norm = matplotlib.colors.Normalize(vmin=np.min(GP),vmax=np.max(GP))
        c_m = matplotlib.cm.Spectral
        s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
        s_m.set_array([])
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')  
        ax.plot_surface(XP,YP,ZP,linewidth=0, rstride=5, cstride=5, facecolors=s_m.to_rgba(GP))
        plt.colorbar(s_m)
        plt.show() 

    """
    # 1) Calculate moving rms and average, see moving()
    # 2) ...
    # 3) ... 
    
    import numpy as np

    if np.asarray(TKO) is not TKO:
        TKO = np.asarray(TKO)
    if np.asarray(AZM) is not AZM:
        AZM = np.asarray(AZM)

    strike = strike * np.pi/180.0
    dip = dip * np.pi / 180.0
    rake = rake * np.pi / 180.0
    gamma = gamma * np.pi / 180.0
    TKO = TKO * np.pi / 180.0
    AZM = AZM * np.pi / 180.0

    Gp = np.cos(TKO)*(np.cos(TKO)*(np.sin(gamma)*(2*np.cos(dip)**2 - (2*sigma)/(2*sigma - 1)) + np.sin(2*dip)*np.cos(gamma)*np.sin(rake)) - np.cos(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) + np.sin(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma))) + np.sin(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.cos(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) + np.sin(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) - np.sin(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) + np.cos(AZM)*np.sin(TKO)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
    Gs = ((np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(gamma)*(2*np.cos(dip)**2 - (2*sigma)/(2*sigma - 1)) + np.sin(2*dip)*np.cos(gamma)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2))))**2 + (np.cos(TKO)*(np.cos(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.sin(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.cos(AZM)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma))))**2)**(1/2)
    Gsh = np.cos(TKO)*(np.cos(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.sin(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.cos(AZM)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)))
    Gsv = np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma)) + np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(gamma)*(2*np.cos(dip)**2 - (2*sigma)/(2*sigma - 1)) + np.sin(2*dip)*np.cos(gamma)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(gamma))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(gamma)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(gamma)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(gamma)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(gamma)) - np.cos(AZM)*np.cos(TKO)*(np.cos(gamma)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(gamma)*((2*sigma)/(2*sigma - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
    
    return (Gp, Gs, Gsh, Gsv)
    
