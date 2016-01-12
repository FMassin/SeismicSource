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

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import obspy.imaging.scripts.mopad





def sheartensile_radpat(mt, TKO, AZM, wave='None', poisson=0.25) : 
    """
    rpgen(strike,dip,rake,isotropic,poisson, TKO, AZM) calculates P-wave, S-wave,
    SH-wave and SV-wave radiation pattern unp.sing shear-tensile source model
    presented in [see references 2, 3, 4 for details]. All input angles 
    (strike, dip, rake of the fault, tensile angle isotropic, takeoff angle 
    TKO and azimuth from the source to the observation point AZM) should 
    be in radians. The takeoff angle is measure from bottom. The function 
    returns matrices of the same size as input TKO and AZM matrices. 

    :type strike, dip, rake: ObsPy : variables
    :param strike, dip, rake: fault plane parameters (radians).
    :type isotropic: variable
    :param isotropic: tensile angle in radians (0 radians for pure shear faulting, 
        pi/2 radians for pure tensile opening).
    :type poisson: variable
    :param poisson:  Poisson's ratio.
    :type TKO: variable
    :param TKO: matrix of takeoff angles for which to calculate the 
        corresponding radiation pattern coefficients (radians, the takeoff 
        angles are measured from bottom).
    :type AZM: variable
    :param AZM: matrix of corresponding azimuths (in radians) for which the 
        radiation pattern coefficients should be calculated.
    :type wave: String
    :param wave: type of wave to compute
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
        
        AZIMUTH = np.arange(0,360,5)*np.pi/180
        TAKEOFF = np.arange(0,180,5)*np.pi/180
        [AZIMUTH,TAKEOFF] = np.meshgrid(AZIMUTH,TAKEOFF)

        strike = 50 * np.pi/180.0
        dip = 60 * np.pi/180.0
        rake = -90 * np.pi/180.0
        poisson = 0.25
        isotropic = 0 * np.pi/180.0
        GP = source.sheartensile_radpat(strike,dip,rake,isotropic,poisson, TAKEOFF, AZIMUTH, 'P')
               

    .. plot::

        strike = 50 * np.pi/180.0
        dip = 60 * np.pi/180.0
        rake = -90 * np.pi/180.0
        poisson = 0.25
        isotropic = 0 * np.pi/180.0
        source.sheartensile_plot(strike,dip,rake,isotropic,poisson, 'P') # try 'P' or 'S' or 'SH' or 'SV' or any ratio

    """
    # 1) Calculate moving rms and average, see moving()
    # 2) ...
    # 3) ... 
    
    if len(mt) == 3 :
        strike, dip, rake = mt
        isotropic = 0
        deviatoric = 0

    elif len(mt) == 4 :
        strike, dip, rake, isotropic = mt
        deviatoric = 0

    elif len(mt) == 5 :
        strike, dip, rake, isotropic, deviatoric = mt

    else:
        print 'Can t yet compute this moment tensor.'

    if np.asarray(TKO) is not TKO:
        TKO = np.asarray(TKO)
    if np.asarray(AZM) is not AZM:
        AZM = np.asarray(AZM)

    if wave in ('None', 'P', 'P-wave', 'P wave'):
        G = np.cos(TKO)*(np.cos(TKO)*(np.sin(isotropic)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(isotropic)*np.sin(rake)) - np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) + np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic))) + np.sin(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) + np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) - np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) + np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
    
    elif wave in ('S', 'S-wave', 'S wave'):
        G = np.sqrt( sheartensile_radpat(mt, TKO, AZM, 'SH', poisson)**2 + sheartensile_radpat(mt, TKO, AZM, 'SV', poisson)**2 )

    elif wave in ('SH', 'Sh', 'Sh-wave', 'Sh wave', 'SH-wave', 'SH wave'):
        G = np.cos(TKO)*(np.cos(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.cos(AZM)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)))
        
    elif wave in ('SV', 'Sv', 'Sv-wave', 'Sv wave', 'SV-wave', 'SV wave'):
        G = np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(isotropic)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(isotropic)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
    
    elif wave in ('S/P', 's/p'):
        G = sheartensile_radpat(mt, TKO, AZM, 'S', poisson)/sheartensile_radpat(mt, TKO, AZM, 'P', poisson)  

    elif wave in ('P/S', 'p/s'):
        G = sheartensile_radpat(mt, TKO, AZM, 'P', poisson)/sheartensile_radpat(mt, TKO, AZM, 'S', poisson)  

    elif wave in ('SH/P', 'sh/p'):
        G = sheartensile_radpat(mt, TKO, AZM, 'SH', poisson)/sheartensile_radpat(mt, TKO, AZM, 'P', poisson)

    elif wave in ('SV/P', 'sv/p'):
        G = sheartensile_radpat(mt, TKO, AZM, 'SV', poisson)/sheartensile_radpat(mt, TKO, AZM, 'P', poisson)

    elif wave in ('SH/S', 'sh/s'):
        G = sheartensile_radpat(mt, TKO, AZM, 'SH', poisson)/sheartensile_radpat(mt, TKO, AZM, 'S', poisson)

    elif wave in ('SV/S', 'sv/s'):
        G = sheartensile_radpat(mt, TKO, AZM, 'SV', poisson)/sheartensile_radpat(mt, TKO, AZM, 'S', poisson)
    else:
        print 'Can t yet compute this wave type.'

    return G


def sheartensile_plot(mt, wave='None', poisson=0.25) :     
        
    AZIMUTH = np.arange(0,360,5)*np.pi/180
    TAKEOFF = np.arange(0,180,5)*np.pi/180
    [AZIMUTH,TAKEOFF] = np.meshgrid(AZIMUTH,TAKEOFF)

    G = sheartensile_radpat(mt, TAKEOFF, AZIMUTH, wave, poisson)

    G[G < -10] = -np.inf
    G[G > 10 ] = np.inf

    scale = np.abs(G)
    Y = scale*np.cos(AZIMUTH) * np.sin(TAKEOFF)
    X = scale*np.sin(AZIMUTH) * np.sin(TAKEOFF)
    Z = scale*-np.cos(TAKEOFF)

    # initializing the colormap machinery
    norm = matplotlib.colors.Normalize(vmin=np.min(G),vmax=np.max(G))
    c_m = matplotlib.cm.Spectral
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')  
    ax.plot_surface(X,Y,Z,linewidth=0, rstride=2, cstride=2, facecolors=s_m.to_rgba(G))
    plt.colorbar(s_m)
    plt.show() 

        
