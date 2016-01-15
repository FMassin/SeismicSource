# -*- coding: utf-8 -*-
"""
source - Module to model seismic sources
======================================================
This module ...

.. rubric:: Example
    
    

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
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import obspy.imaging.scripts.mopad


def sphere(r=1.,n=100.):
    """
    produce the polar coordinates of a sphere. 

    :type r, n: variables
    :param r, n: radius and number of resolution points.
    
    :rtype : list
    :return : 3d list of phi, theta and radius.

    .. seealso::

        numpy.linspace : produces the resolution vectors.

        numpy.meshgrid : produce the grid from the vectors.

    .. rubric:: Example

        # 100 coordinates on sphere or radius 5
        points = source.sphere_spherical(r=5, n=100)
        # Unit sphere of 50 points
        points = source.sphere_spherical(n=50)               

    .. plot::

        # run one of the example
        points = spherical_to_cartesian(points)
        fig = plt.figure()
        ax=fig.gca(projection='3d')
        ax.set_aspect("equal")
        ax.plot_wireframe(points[0], points[1], points[2], color="r")
        plt.show() 

    """
    # Get the resolution (sorry this is ugly)
    c = 0.038 ; 
    na = n**(.5+c) 
    nt = n**(.5-c)

    [AZIMUTH,TAKEOFF] = np.meshgrid( np.linspace(0, 2*np.pi, na), np.linspace(0, np.pi, nt) )

    RADIUS = np.ones(AZIMUTH.shape)*r

    return [ AZIMUTH, TAKEOFF, RADIUS ]


def cartesian_to_spherical(vector):
    """This file is part of the program relax (Edward d'Auvergne).
    
    http://svn.gna.org/svn/relax/1.3/maths_fns/coord_transform.py

    Convert the Cartesian vector [x, y, z] to spherical coordinates [theta, phi, r].

    The parameter r is the radial distance, theta is the polar angle, and phi is the azimuth.


    @param vector:  The Cartesian vector [x, y, z].
    @type vector:   numpy rank-1, 3D array
    @return:        The spherical coordinate vector [theta, phi, r].
    @rtype:         numpy rank-1, 3D array
    """

    # Make sure we got np.array 
    if np.asarray(vector) is not vector:
        vector = np.asarray(vector)

    # The radial distance.
    r = np.sqrt((vector**2).sum(axis=0))

    # Unit vector.
    unit = vector / r

    # The polar angle.
    theta = np.arccos(unit[2])

    # The azimuth.
    phi = np.arctan2(unit[1], unit[0])

    # Return the spherical coordinate vector.
    return [theta, phi, r]


def spherical_to_cartesian(vector):
    """This file is part of the program relax (Edward d'Auvergne).
    
    http://svn.gna.org/svn/relax/1.3/maths_fns/coord_transform.py

    Convert the spherical coordinate vector [theta, phi, r] to the Cartesian vector [x, y, z].

    The parameter r is the radial distance, theta is the polar angle, and phi is the azimuth.


    @param vector:  The spherical coordinate vector [phi, theta, r].
    @type vector:   3D array or list
    @param cart_vect:       The Cartesian vector [x, y, z].
    @type cart_vect:        3D array or list
    """

    # Trig alias.
    sin_theta = np.sin(vector[1])

    # Unit vector if r is missing
    if len(vector) == 2 :
        r=1
    else:
        r=vector[2]

    # The vector.
    x = r * np.cos(vector[0]) * sin_theta
    y = r * np.sin(vector[0]) * sin_theta
    z = r * np.cos(vector[1])

    return [x, y, z]


def fullmt(mt):
    """takes 6 comp moment tensor and returns full 3x3 moment tensor"""
    mt_full = np.array(([[mt[0], mt[3], mt[4]],
                         [mt[3], mt[1], mt[5]],
                         [mt[4], mt[5], mt[2]]]))
    return mt_full



class Model(object):

    def __init__(self, mt, poisson=0.25):
        self.mt = mt
        self.poisson = poisson

    def farfield_radpat(self, wave='P', obs_cart='None', obs_sph=sphere(r=1.,n=500.)):
        """
        Returns the P farfield radiation pattern
        based on Aki & Richards Eq 4.29

        :param mt: Focal mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - the
                   six independent components of the moment tensor)

        :type wave: String
        :param wave: type of wave to compute

        :param points: 3D vector array with shape [3,npts] (x,y,z) or [2,npts]
                       (theta,phi) The normalized displacement of the moment
                       tensor source is computed at these points.

        :return: 3D vector array with shape [3,npts] that contains the
                 displacement vector for each grid point
        """
        
        # Get full mt
        ## ... from 6-elt mt
        Mpq = fullmt(self.mt)
        ## ... from strike, dip, rake, iso, clvd
        ## ... from strike, dip, rake

        # Get observation points
        ## Get unit sphere, or spherical coordinate if given 
        if obs_cart == 'None' :
            obs_cart = spherical_to_cartesian(obs_sph)
        ## Make sure they are np.array 
        if np.asarray(obs_cart) is not obs_cart:
            obs_cart = np.asarray(obs_cart)
        ## Keeping that in mind
        requestdimension = obs_cart.shape
        obs_cart = np.reshape(obs_cart, (3, np.prod(requestdimension)/3))
        
        # Displacement array    
        # precompute directional cosine array
        dists = np.sqrt(obs_cart[0] * obs_cart[0] + obs_cart[1] * obs_cart[1] +
                        obs_cart[2] * obs_cart[2])

        # In gamma, all points are taken to a unit distance, same angle:
        gammas = obs_cart / dists

        # initialize displacement array
        ndim, npoints = obs_cart.shape
        disp = np.empty(obs_cart.shape)

        # loop through points
        for ipoint in range(npoints):
            
            gamma = gammas[:, ipoint]

            if wave in ('S', 'S wave', 'S-wave'):

                Mp = np.dot(Mpq, gamma)

                # loop through displacement component [n index]
                for n in range(ndim):
                    psum = 0.0
                    for p in range(ndim):
                        deltanp = int(n == p)
                        psum += (gamma[n] * gamma[p] - deltanp) * Mp[p]
                    disp[n, ipoint] = psum

            elif wave in ('P', 'P wave', 'P-wave'):

                gammapq = np.outer(gamma, gamma)
                gammatimesmt = gammapq * Mpq
                
                # loop through displacement component [n index]
                for n in range(ndim):
                    disp[n, ipoint] = gamma[n] * np.sum(gammatimesmt.flatten())

        # Reshape to request dimensions
        obs_cart = np.reshape(obs_cart, requestdimension)
        disp = np.reshape(disp, requestdimension)

        return disp, obs_cart

    def farfield_surf(self, wave='P',style='None') :                     
        

        # Get radiation pattern
        G, pointscart = self.farfield_radpat(wave)

        AZIMUTH, TAKEOFF, RADIUS = cartesian_to_spherical(pointscart)

        # Take a cartesian radpat, make a spherical (much easier to plot)
        magn = np.sum(G * pointscart, axis=0)
        magn /= np.max(np.abs(magn))
        G = np.sqrt( G[0]**2 + G[1]**2 + G[2]**2 ) * magn
        G = np.reshape(G, AZIMUTH.shape )

        # Styles
        if style in ('sign', 'bb', 'beachball', 'DC', 'fp'):
            G[G < 0] = -1
            G[G >= 0] = 1
            scale=1
        elif style == 'None' : 
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
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show() 





 



def radpat_sheartensile_vavryeuk(mt, wave='P', TKO=np.arange(0,181,1)*np.pi/180, AZM=np.arange(0,361,1)*np.pi/180, poisson=0.25) : 
    """
    I don't understand how the rake is used in this function
    It seems to be a little different than farfield_radpat 

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
    
    # Moment tensor
    ## Getting various formats
    if len(mt) == 3 :
        strike, dip, rake = mt
        isotropic = 0
        deviatoric = 0

    elif len(mt) == 4 :
        strike, dip, rake, isotropic = mt
        deviatoric = 0

    elif len(mt) == 5 :
        strike, dip, rake, isotropic, deviatoric = mt

    elif len(mt) == 6 :
        print 'We need to find a piece of code to make that...'

    elif len(mt) == 9 :
        print 'We need to find a piece of code to make that...'

    else:
        print 'Can t yet compute this moment tensor.'

    ## convert isotropic ratio to angle
    isotropic = np.arcsin(isotropic)
    
    # Spherical angles
    ## Make sure we got np.array 
    if np.asarray(TKO) is not TKO:
        TKO = np.asarray(TKO)
    if np.asarray(AZM) is not AZM:
        AZM = np.asarray(AZM)

    ## Make a 2D grid if 2 1d arrays
    if np.shape(AZM.shape) == np.shape(TKO.shape) and np.shape(AZM.shape) == (1,) and AZM.shape != TKO.shape:
        [AZM,TKO] = np.meshgrid(AZM,TKO)


    # Radiation patterns
    ## Tensile definitions by Vavryèuk (2001)
    if wave in ('P', 'P-wave', 'P wave'):
        G = np.cos(TKO)*(np.cos(TKO)*(np.sin(isotropic)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(isotropic)*np.sin(rake)) - np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) + np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic))) + np.sin(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) + np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) - np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) + np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
    
    elif wave in ('SH', 'Sh', 'Sh-wave', 'Sh wave', 'SH-wave', 'SH wave'):
        G = np.cos(TKO)*(np.cos(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.cos(AZM)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)))
        
    elif wave in ('SV', 'Sv', 'Sv-wave', 'Sv wave', 'SV-wave', 'SV wave'):
        G = np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(isotropic)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(isotropic)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
    
    ## Re-using the same programme to get other things ... 
    elif wave in ('S', 'S-wave', 'S wave'):
        G = np.sqrt( radpat_sheartensile_vavryeuk(mt, 'SH', TKO, AZM, poisson)**2 + radpat_sheartensile_vavryeuk(mt, 'SV', TKO, AZM, poisson)**2 )

    elif wave in ('S/P', 's/p'):
        G = radpat_sheartensile_vavryeuk(mt, 'S', TKO, AZM, poisson)/radpat_sheartensile_vavryeuk(mt, 'P', TKO, AZM, poisson)  

    elif wave in ('P/S', 'p/s'):
        G = radpat_sheartensile_vavryeuk(mt, 'P', TKO, AZM, poisson)/radpat_sheartensile_vavryeuk(mt, 'S', TKO, AZM, poisson)  

    elif wave in ('SH/P', 'sh/p'):
        G = radpat_sheartensile_vavryeuk(mt, 'SH', TKO, AZM, poisson)/radpat_sheartensile_vavryeuk(mt, 'P', TKO, AZM, poisson)

    elif wave in ('SV/P', 'sv/p'):
        G = radpat_sheartensile_vavryeuk(mt, 'SV', TKO, AZM, poisson)/radpat_sheartensile_vavryeuk(mt, 'P', TKO, AZM, poisson)

    elif wave in ('SH/S', 'sh/s'):
        G = radpat_sheartensile_vavryeuk(mt, 'SH', TKO, AZM, poisson)/radpat_sheartensile_vavryeuk(mt, 'S', TKO, AZM, poisson)

    elif wave in ('SV/S', 'sv/s'):
        G = radpat_sheartensile_vavryeuk(mt, 'SV', TKO, AZM, poisson)/radpat_sheartensile_vavryeuk(mt, 'S', TKO, AZM, poisson)

    ## Making sure you get that error.
    else:
        print 'Can t yet compute this wave type.'

    return G


def surfacemotions(mt, wave='P', poisson=0.25, style='None') :     
        
    AZIMUTH = np.arange(0,365,5)*np.pi/180
    TAKEOFF = np.arange(0,185,5)*np.pi/180
    [AZIMUTH,TAKEOFF] = np.meshgrid(AZIMUTH,TAKEOFF)


    G = radpat_sheartensile_vavryeuk(mt, wave, TAKEOFF, AZIMUTH, poisson)
    
    #G[G < -10] = -np.inf
    #G[G > 10 ] = np.inf

    if style in ('sign', 'bb', 'beachball', 'DC', 'fp'):
        G[G < 0] = -1
        G[G >= 0] = 1
        scale=1
    elif style == 'None' : 
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
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show() 

def isotropiceffectonSP(mt=[0, 0, 0, 0]):

    rms=[]
    pythagoresum=[]
    average=[]

    for x in range(0, 100):
        mt[3] = x / 100.
        Prd = radpat_sheartensile_vavryeuk(mt,'P') 
        Srd = radpat_sheartensile_vavryeuk(mt,'S') 
        P = np.sum(Prd**2)/np.prod(Prd.shape) # classic rms
        S = np.sum(Srd**2)/np.prod(Srd.shape)
        rms.append(S/P)
        P = np.sqrt(np.sum(Prd**2)) # pythagorian sum
        S = np.sqrt(np.sum(Srd**2)) 
        pythagoresum.append(S/P)
        P = np.average(np.abs(Prd)) # amplitude average
        S = np.average(np.abs(Srd)) 
        average.append(S/P)        

    plt.plot(np.arange(100),rms, 'r--', np.arange(100),pythagoresum, 'bs', np.arange(100),average, 'g^')
    plt.xlabel('Isotropic %')
    plt.ylabel('S/P')
    plt.title('This how S/P ratio evolve with isotropic component')
    plt.grid(True)
    plt.legend( ["RMS", "Pythagorian sum", "Average" ] )
    plt.show()






class SeismicSource(object):

    # general attribut definition
    notes = 'Ceci est a moi, A MOUHA, A-MOU-HAAAAAA'

    def __init__(self, mt):

        self.mt = mt       

    def surfacemotions(self, wave='P', poisson=0.25, style='None'):

        surfacemotions(self.mt, wave, poisson, style) 

    def isotropiceffectonSP(self):

        isotropiceffectonSP(self.mt)

        
