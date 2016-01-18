# -*- coding: utf-8 -*-
"""
source - Module for seismic sources modeling
======================================================
This module compute the seismic wave displacement 
using the theoritical definitions of the seismic 
source.

.. rubric:: Example
    
.. figure:: /_images/test.png

.. note::

    For ...

:copyright:
    The ...
:license:
    ...
"""


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import obspy.imaging.scripts.mopad


def mt_full(mt):
    """takes 6 comp moment tensor and returns full 3x3 moment tensor"""
    mt = np.array(([[mt[0], mt[3], mt[4]],
        [mt[3], mt[1], mt[5]],
        [mt[4], mt[5], mt[2]]]))
    return mt


def mt_angles(mt): 

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
        strike, dip, rake, isotropic, deviatoric, poubelle = mt
        print 'We need to find a piece of code to make that...'

    elif len(mt) == 9 :
        print 'We need to find a piece of code to make that...'

    else:
        print 'Can t yet compute this moment tensor.'

    return [strike, dip, rake, isotropic, deviatoric]


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
    return [phi, theta, r]


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
    y = r * np.cos(vector[0]) * sin_theta
    x = r * np.sin(vector[0]) * sin_theta
    z = r * -np.cos(vector[1])

    return [x, y, z]


def plot_seismicsourcemodel(G, XYZ, style='None') : 

    magn = np.sum(G * XYZ, axis=0)
    magn /= np.max(np.abs(magn))

    amplitudes = np.sqrt( G[0]**2 + G[1]**2 + G[2]**2 ) * magn

    # Styles
    if style in ('sign', 'bb', 'beachball', 'DC', 'fp'):
        G[G < 0] = -1
        G[G >= 0] = 1
        scale=1
    elif style == 'None' : 
        scale = np.abs(amplitudes)
    
    XYZ[0] *= scale 
    XYZ[1] *= scale 
    XYZ[2] *= scale

    # initializing the colormap machinery
    norm = matplotlib.colors.Normalize(vmin=np.min(amplitudes),vmax=np.max(amplitudes))
    c_m = matplotlib.cm.Spectral
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')  
    ax.plot_surface(XYZ[0], XYZ[1], XYZ[2],linewidth=0, rstride=1, cstride=1, facecolors=s_m.to_rgba(amplitudes), alpha=0.5)
    plt.colorbar(s_m)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show() 


def energy_seismicsourcemodel(G, XYZ) :    

        #amplitudes_correction = np.sum(G * XYZ, axis=0)
        #amplitudes_correction /= np.max(np.abs(amplitudes_correction))
        amplitudes = np.sqrt( G[0]**2 + G[1]**2 + G[2]**2 ) #* amplitudes_correction            
        
        # Classic rms
        rms = np.sum(amplitudes**2)/np.prod(amplitudes.shape)  
        # Euclidian norm
        norm = np.sqrt(np.sum(amplitudes**2))  
        # Amplitude average
        average = np.average(np.abs(amplitudes))  

        return [rms, norm, average]


class Aki_Richards(object):

    def __init__(self, mt, poisson=0.25):
        self.mt = mt
        self.poisson = poisson
        
    def radpat(self, wave='P', obs_cart='None', obs_sph=sphere(r=1.,n=1000.)):
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
        Mpq = mt_full(self.mt)

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

    def plot(self, wave='P',style='None') :    

        # Get radiation pattern
        G, XYZ = self.radpat(wave)
        plot_seismicsourcemodel(G, XYZ, style='None')

    def energy(self, wave='P') :    

        # Get radiation pattern and estimate energy
        G, XYZ = self.radpat(wave)  
        estimators = energy_seismicsourcemodel(G, XYZ)
        return estimators


class Vavryeuk(object):

    def __init__(self, mt, poisson=0.25):
        self.mt = mt
        self.poisson = poisson
        
    def radpat(self, wave='P', obs_cart=spherical_to_cartesian(sphere(r=1.,n=1000.)), obs_sph='None'):
        #(mt, wave='P', TKO=np.arange(0,181,1)*np.pi/180, AZM=np.arange(0,361,1)*np.pi/180, poisson=0.25) : 
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
        
        poisson = self.poisson

        # Moment tensor
        strike, dip, rake, isotropic, deviatoric = mt_angles(self.mt)        

        ## convert isotropic ratio to angle
        isotropic = np.arcsin(isotropic)
        

        # Get observation points
        ## Get unit sphere, or spherical coordinate if given 
        if obs_sph == 'None' :
            obs_sph = cartesian_to_spherical(obs_cart)
        else:
            obs_cart = spherical_to_cartesian(obs_sph)

        ## Make sure they are np.array 
        if np.asarray(obs_sph) is not obs_sph:
            obs_sph = np.asarray(obs_sph)
        if np.asarray(obs_cart) is not obs_cart:
            obs_cart = np.asarray(obs_cart)

        
        # Radiation patterns
        # G is the amplitude for each observation point
        ## Observations are given in spherical angles
        AZM = obs_sph[0]
        TKO = obs_sph[1]

        ## Make sure we got np.array 
        if np.asarray(TKO) is not TKO:
            TKO = np.asarray(TKO)
        if np.asarray(AZM) is not AZM:
            AZM = np.asarray(AZM)

        ## Tensile definitions by Vavryèuk (2001)
        if wave in ('P', 'P-wave', 'P wave'):
            G = np.cos(TKO)*(np.cos(TKO)*(np.sin(isotropic)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(isotropic)*np.sin(rake)) - np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) + np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic))) + np.sin(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) + np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) - np.sin(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) + np.cos(AZM)*np.sin(TKO)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
        
        elif wave in ('SH', 'Sh', 'Sh-wave', 'Sh wave', 'SH-wave', 'SH wave'):
            G = np.cos(TKO)*(np.cos(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.cos(AZM)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)))
            
        elif wave in ('SV', 'Sv', 'Sv-wave', 'Sv wave', 'SV-wave', 'SV wave'):
            G = np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(isotropic)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(isotropic)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
        
        ## Re-using the same programme to get other things ... 
        elif wave in ('S', 'S-wave', 'S wave'):

            # for such definition this the less ugly
            Gsh = np.cos(TKO)*(np.cos(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.cos(AZM)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)))
            Gsv = np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic)) + np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(isotropic)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(isotropic)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(isotropic))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(isotropic)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(isotropic)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(isotropic)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(isotropic)) - np.cos(AZM)*np.cos(TKO)*(np.cos(isotropic)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(isotropic)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))

            G = np.sqrt(Gsh**2 + Gsv**2)

        elif wave in ('S/P', 's/p'):
            G = self.radpat('S', )/self.radpat('P', obs_sph = obs_sph)  

        elif wave in ('P/S', 'p/s'):
            G, poubelle  = self.radpat('P', obs_sph = obs_sph)/self.radpat('S', obs_sph = obs_sph)  

        elif wave in ('SH/P', 'sh/p'):
            G, poubelle  = self.radpat('SH', obs_sph = obs_sph)/self.radpat('P', obs_sph = obs_sph)

        elif wave in ('SV/P', 'sv/p'):
            G, poubelle  = self.radpat('SV', obs_sph = obs_sph)/self.radpat('P', obs_sph = obs_sph)

        elif wave in ('SH/S', 'sh/s'):
            G, poubelle  = self.radpat('SH', obs_sph = obs_sph)/self.radpat('S', obs_sph = obs_sph)

        elif wave in ('SV/S', 'sv/s'):
            G, poubelle  = self.radpat('SV', obs_sph = obs_sph)/self.radpat('S', obs_sph = obs_sph)

        ## Making sure you get that error.
        else:
            print 'Can t yet compute this wave type.'

        ## transform G into vector x,y,z          
        G_cart = spherical_to_cartesian(np.asarray([AZM, TKO, G]))
        obs_cart = spherical_to_cartesian(np.asarray([AZM, TKO]))

        #return G, [AZM, TKO] 
        return np.asarray(G_cart), np.asarray(obs_cart)


    def plot(self, wave='P',style='None') :    

        # Get radiation pattern and plot
        G, XYZ = self.radpat(wave)
        plot_seismicsourcemodel(G, XYZ, style='None')

    def energy(self, wave='P') :    

        # Get radiation pattern and estimate energy
        G, XYZ = self.radpat(wave)  
        estimators = energy_seismicsourcemodel(G, XYZ)
        return estimators



class SeismicSource(object):

    # general attribut definition
    notes = 'Ceci est a moi, A MOUHA, A-MOU-HAAAAAA'

    def __init__(self, mt, poisson=0.25):

        self.Aki_Richards = Aki_Richards(mt,poisson)  
        self.Vavryeuk = Vavryeuk(mt,poisson)  


# function skeleton
# def ggg(...):
#         """
#         Run a ...
#
#         :param type: String that specifies ... (e.g. ``'skeleton'``).
#         :param options: Necessary keyword arguments for the respective
#             type.
#
#         .. note::
#
#             The raw data is/isn't accessible anymore afterwards.
#
#         .. rubric:: _`Supported type`
#
#         ``'skeleton'``
#             Computes the ... (uses 
#               :func:`dependency.core.Class1.meth1`).
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





        
