# -*- coding: utf-8 -*-
"""
source - Module for seismic sources modeling
______________________________________________________________________
This module provide class hierarchy for earthquake modeling and 
representation.

The seismic wave displacement  is computed using the theoritical 
definitions of the seismic source.
"""


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from obspy.imaging.scripts.mopad import MomentTensor

# Functions and classes are ordered from general to specific

def sphere(r=1.,n=100.):
    """
    Produce the polar coordinates of a sphere. 
    ______________________________________________________________________
    :type r, n: variables
    :param r, n: radius and number of resolution points.
    
    :rtype : list
    :return : 3d list of azimuth, polar angle and radial distance.

    .. seealso::

        numpy.linspace : produces the resolution vectors.

        numpy.meshgrid : produce the grid from the vectors.

    .. rubric:: Example

        # import the module
        import source

        # 100 coordinates on sphere or radius 5
        points = source.sphere(r=5, n=100)

        # Unit sphere of 50 points
        points = source.sphere(n=50)               

    .. plot::

        # run one of the example
        points = source.spherical_to_cartesian(points)

        # plot using matplotlib 
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax=fig.gca(projection='3d')
        ax.set_aspect("equal")
        ax.plot_wireframe(points[0], points[1], points[2], color="r")
        plt.show() 
    ______________________________________________________________________
    """
    # Get the resolution (sorry this is ugly)
    c = 0.038 ; 
    na = n**(.5+c) 
    nt = n**(.5-c)

    [AZIMUTH,TAKEOFF] = np.meshgrid( np.linspace(0, 2*np.pi, na+1), np.linspace(0, 1*np.pi, nt) )

    RADIUS = np.ones(AZIMUTH.shape)*r

    return [ AZIMUTH, TAKEOFF, RADIUS ]


def cartesian_to_spherical(vector):
    """
    Convert the Cartesian vector [x, y, z] to spherical coordinates 
    [azimuth, polar angle, radial distance].
    ______________________________________________________________________
    :type vector : 3D array, list | np.array
    :param vector :  The vector of cartessian coordinates.

    :rtype : 3D array, np.array
    :return : The spherical coordinate vector.
    
    .. note::

        This file is extracted & modified from the program relax (Edward 
            d'Auvergne).

    .. seealso::
    
        http://svn.gna.org/svn/relax/1.3/maths_fns/coord_transform.py
    ______________________________________________________________________
    
    """

    # Make sure we got np.array 
    if np.asarray(vector) is not vector:
        vector = np.asarray(vector)

    # The radial distance.
    radius = np.sqrt((vector**2).sum(axis=0))

    # The horizontal radial distance.
    rh = np.sqrt(vector[0]**2 + vector[1]**2) 

    # The polar angle.
    takeoff = np.arccos( vector[2] / radius )
    takeoff[radius == 0.0] = np.pi / 2 * np.sign(vector[2][radius == 0.0])
    #theta = np.arctan2(vector[2], rh)

    # The azimuth.
    azimuth_trig = np.arctan2(vector[1], vector[0])

    # Return the spherical coordinate vector.
    return [azimuth_trig, takeoff, radius]


def spherical_to_cartesian(vector):
    """
    Convert the spherical coordinates [azimuth, polar angle
    radial distance] to Cartesian coordinates [x, y, z].

    ______________________________________________________________________
    :type vector : 3D array, list | np.array
    :param vector :  The spherical coordinate vector.

    :rtype : 3D array, np.array
    :return : The vector of cartesian coordinates.
    
    .. note::

        This file is extracted & modified from the program relax (Edward 
            d'Auvergne).

    .. seealso::
    
        http://svn.gna.org/svn/relax/1.3/maths_fns/coord_transform.py
    ______________________________________________________________________
    """

    # Unit vector if r is missing
    if len(vector) == 2 :
        radius =1
    else:
        radius=vector[2]

    # Trig alias.
    sin_takeoff = np.sin(vector[1])

    # The vector.
    x = radius * sin_takeoff * np.cos(vector[0])
    y = radius * sin_takeoff * np.sin(vector[0])
    z = radius * np.cos(vector[1])


    return [x, y, z]


def project_vectors(b, a):
    """
    Project the vectors b on vectors a.
    ______________________________________________________________________
    :type b : 3D np.array
    :param b :  The cartesian coordinates of vectors to be projected.

    :type a : 3D np.array
    :param a :  The cartesian coordinates of vectors to project onto.

    :rtype : 3D np.array
    :return : The cartesian coordinates of b projected on a.
    
    .. note::

        There is maybe better (built in numpy) ways to do this.

    .. seealso::
    
        :func:`numpy.dot()` 
        :func:`numpy.dotv()`
    ______________________________________________________________________
    """

    # Project 
    ## ratio (norm of b on a / norm of a)
    proj_norm = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] 
    proj_norm /=  (a[0]**2 + a[1]**2 + a[2]**2)+0.0000000001
    ## project on a
    b_on_a = proj_norm * a


    return b_on_a


def vector_normal(XYZ, v_or_h):
    """
    Compute the vertical or horizontal normal vectors.
    ______________________________________________________________________
    :type XYZ : 3D np.array
    :param b :  The cartesian coordinates of vectors.

    :type v_or_h : string
    :param v_or_h :  The cartesian coordinates of vectors.

    :rtype : 3D np.array
    :return : Requested normal on vertical or horizontal plane 
        ('v|vert|vertical' or 'h|horiz|horizontal').
    
    .. note::

        There is maybe better (built in numpy) ways to do this.

    .. seealso::
    
        :func:`numpy.cross()` 
    ______________________________________________________________________
    """

    oneheightydeg =  np.ones(XYZ.shape) * np.pi
    ninetydeg =  np.ones(XYZ.shape) * np.pi/2.

    ATR = np.asarray(cartesian_to_spherical(XYZ))

    if v_or_h in ('n', 'nhr', 'normal horizontal radial', 'norm horiz rad'): 
        ATR_n = np.array([ ATR[0] , ATR[1] - ninetydeg[0], ATR[2]])
        XYZ_n = np.asarray(spherical_to_cartesian(ATR_n))
    elif v_or_h in ('h', 'horizontal', 'horiz'): 
        ATR_n = np.array([ ATR[0] - ninetydeg[0], ATR[1]+(ninetydeg[0]-ATR[1]) , ATR[2]])
        XYZ_n = np.asarray(spherical_to_cartesian(ATR_n))
    elif v_or_h in ('v', 'vertical', 'vertical'): 
        ATR_n = np.array([ ATR[0], ATR[1]-ATR[1] , ATR[2]])
        XYZ_n = np.asarray(spherical_to_cartesian(ATR_n))
    elif v_or_h in ('r', 'radial', 'self'): 
        XYZ_n = XYZ

    return XYZ_n

def mt_full(mt):
    """    
    Takes 6 components moment tensor and returns full 3x3 moment 
    tensor.
    ______________________________________________________________________
    :type mt : list or np.array.
    :param mt : moment tensor NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz, 
        the six independent components).

    :rtype : np.array 3x3.
    :return : full 3x3 moment tensor.

    .. note::

        Use obspy.imaging.scripts.mopad to get correct input.
    ______________________________________________________________________
    """

    # Make sure we got np.array 
    if np.asarray(mt) is not mt:
        mt = np.asarray(mt)

    if len(mt) == 6:        
        mt = np.array(([[mt[0], mt[3], mt[4]],
            [mt[3], mt[1], mt[5]],
            [mt[4], mt[5], mt[2]]]))

    if mt.shape != (3,3) :       
        raise Exception('I/O dimensions: only 1x6 or 3x3 input supported.')
    
    return mt


def mt_angles(mt): 
    """    
    Takes 6 components and returns fps tri-angles in degrees, with 
    deviatoric, isotropic, DC and CLVD percentages.
    ______________________________________________________________________
    :type mt : list or np.array.
    :param mt : moment tensor NM x 6 (Mxx, Myy, Mzz, 
        Mxy, Mxz, Myz, the six independent 
        components).

    :rtype : np.array [[1x3], [1x4]].
    :return : full 3x3 moment tensor.

    .. note::
        Nan value are returned for unknown parameters.

        The deviatoric percentage is the sum of DC and CLVD percentages.

        Use obspy.imaging.scripts.mopad to get correct input.
    ______________________________________________________________________
    """

    # Make sure we got np.array 
    if np.asarray(mt) is not mt:
        mt = np.asarray(mt)

    # Getting various formats
    ## if given strike, dip, rake
    if mt.shape == (3,):
        strike, dip, rake = mt
        DC = 100
        CLVD = 0
        iso = 0
        devi = 100
    
    ## if given [[strike, dip, rake], [strike, dip, rake]] (e.g. by MoPad)
    elif mt.shape == (2,3) :
        strike, dip, rake = mt[0]
        DC = 100
        CLVD = 0
        iso = 0
        devi = 100
    
    ## if given [strike, dip, rake, devi]
    elif mt.shape == (4,):
        strike, dip, rake, devi = mt
        DC = np.nan
        CLVD = np.nan
        iso = 0

    ## if given [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
    elif mt.shape == (6,) :
        
        mt = MomentTensor(mt,'XYZ') 
        
        DC = mt.get_DC_percentage()
        CLVD = mt.get_CLVD_percentage()
        iso = mt.get_iso_percentage()
        devi = mt.get_devi_percentage()

        mt = mt_angles(mt.get_fps())
        strike, dip, rake = mt[0]

    ## if given full moment tensor
    elif mt.shape == (3,3) :

        mt = mt_angles([mt[0,0], mt[1,1], mt[2,2], mt[0,1], mt[0,2], mt[1,2]])

        strike, dip, rake = mt[0]
        DC, CLVD, iso, devi = mt[1]

    else:         
        raise Exception('I/O dimensions: only [1|2]x3, 1x[4|6] and 3x3 inputs supported.')

    return np.array([[strike, dip, rake], [DC, CLVD, iso, devi]])


def plot_seismicsourcemodel(disp, xyz, style='*', mt=None, comp=None, ax=None) : 
    """
    Plot the given seismic wave radiation pattern as a color-coded surface 
    or focal sphere (not exactly as a beach ball diagram).
    ______________________________________________________________________
    :type G : 3D array, list | np.array
    :param G : The vector of cartessian coordinates of the radiation 
        pattern.

    :type xyz : 3D array, list | np.array
    :param xyz : The cartessian coordinates of the origin points of the
        radiation pattern.

    :type style : string
    :param style : type of plot.

    :type mt : list | np.array
    :param mt : moment tensor definition supported by MoPad, used in title.

    :rtype : graphical object
    :return : The vector of cartessian coordinates.

    .. rubric:: _`Supported style`

        ``'*'``
            Plot the amplitudes as surface coding absolute amplitude as 
            radial distance from origin and amplitude polarity with color.
             (uses :func:`numpy.abs`).

        ``beachball``
            Plot a unit sphere, amplitude sign being coded with color.
             (uses :func:`numpy.abs`).

        ``'wireframe'``
            Plot the amplitudes as a mesh, coding absolute amplitude as 
            radial distance. Amplitude polarity is not represented.

        ``'quiver'``
            Plot each vectors of the radiation pattern, sign and amplidude
            being coded with colors. 

        ``'surface'``
            Plot ...

        ``'polarities'``
            Plot ...
    
    .. note::

        The radition pattern is rendered as is.

    .. seealso::
    
        For earthquake focal mechanism rendering, use 
        obspy.imaging.beachball.Beachball() : 
        https://docs.obspy.org/packages/obspy.imaging.html#beachballs

        For more info see  and obspy.imaging.mopad_wrapper.Beachball().

    .. todo::

        make work in given axe
    ______________________________________________________________________
    """

    # Component
    ## Get direction(s)
    if comp == None: 
        amplitude_direction = disp      
    else :
        amplitude_direction = vector_normal(xyz, comp)
    ## Project 
    disp_projected = project_vectors(disp, amplitude_direction) 
    ## Norms
    amplitudes = np.sqrt(disp_projected[0]**2 + disp_projected[1]**2 + disp_projected[2]**2) 
    ## Amplitudes (norm with sign)
    amplitudes *= np.sign(np.sum(disp * amplitude_direction, axis=0)+0.00001)
    ## Easiness
    U, V, W = disp_projected
    X, Y, Z = xyz
    amplitudes_surf = np.abs(amplitudes)/np.max(np.abs(amplitudes))

    # Initializing 
    ## Initializing the colormap machinery
    norm = matplotlib.colors.Normalize(vmin=np.min(amplitudes),vmax=np.max(amplitudes))
    c_m = matplotlib.cm.Spectral
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    ## Initializing the plot
    if ax == None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')  
    plt.xlabel('X')
    plt.ylabel('Y')
    ## Initializing the colorbar
    if style not in ('frame', 'wireframe') : 
        cbar = plt.colorbar(s_m,orientation="horizontal",fraction=0.07)
        cbar.ax.set_xlabel('Raw amplitudes')
    ## Initializing figure keys
    if mt is not None : 
        [strike, dip, rake], [DC, CLVD, iso, deviatoric] = mt_angles(mt) 
        plt.title('SDR [' + str(int(strike)) + ', ' + str(int(dip)) + ', ' + str(int(rake)) + '] IDC [' + str(int(iso)) + ', ' + str(int(DC)) + ', ' + str(int(CLVD)) +']%' )
    
    ## Force axis equal
    arrow_scale = (np.max(xyz)-np.min(xyz))/10
    if style in ('p', 'polarities','b', 'bb', 'beachball'):
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0
    elif style in ('f', 'w', 'frame', 'wireframe', 's', 'surf', 'surface'): 
        max_range = np.array([(X*np.abs(amplitudes)).max()-(X*np.abs(amplitudes)).min(), (Y*np.abs(amplitudes)).max()-(Y*np.abs(amplitudes)).min(), (Z*np.abs(amplitudes)).max()-(Z*np.abs(amplitudes)).min()]).max() / 2.0
    else : 
        max_range = np.array([(X+arrow_scale).max()-(X+arrow_scale).min(), (Y+arrow_scale).max()-(Y+arrow_scale).min(), (Z+arrow_scale).max()-(Z+arrow_scale).min()]).max() / 2.0
    mean_x = X.mean()
    mean_y = Y.mean()
    mean_z = Z.mean()
    ax.set_xlim(mean_x - max_range, mean_x + max_range)
    ax.set_ylim(mean_y - max_range, mean_y + max_range)
    ax.set_zlim(mean_z - max_range, mean_z + max_range)


    # Simple styles   
    ## For a wireframe surface representing amplitudes
    if style in ('f', 'w', 'frame', 'wireframe') :
        
        ax.plot_wireframe(X*amplitudes_surf, Y*amplitudes_surf, Z*amplitudes_surf, rstride=1, cstride=1, linewidth=0.5, alpha=0.5)

    ## For focal sphere, with amplitude sign (~not a beach ball diagram) on unit sphere 
    if style in ('*', 'p', 'polarities'):

        polarity_area = amplitudes.copy()
        polarity_area[amplitudes > 0] = np.nan  
        polarity_area[amplitudes <= 0] = 1
        ax.plot_wireframe(X*polarity_area, Y*polarity_area, Z*polarity_area, color='r', rstride=1, cstride=1, linewidth=.5, alpha=0.5)

        polarity_area[amplitudes <= 0] = np.nan 
        polarity_area[amplitudes > 0] = 1
        ax.plot_wireframe(X*polarity_area, Y*polarity_area, Z*polarity_area, color='b', rstride=1, cstride=1, linewidth=.5, alpha=0.5)

    ## For ~ beach ball diagram 
    if style in ('b', 'bb', 'beachball'):

        polarity_area = amplitudes.copy()
        polarity_area[amplitudes >= 0] = 1  
        polarity_area[amplitudes < 0] = -1
        ax.plot_surface(X, Y, Z, linewidth=0, rstride=1, cstride=1, facecolors=s_m.to_rgba(polarity_area)) 


    # Complexe styles
    ## For arrow vectors on obs points
    if style in ('*', 'q', 'a', 'v', 'quiver', 'arrow', 'vect', 'vector', 'vectors') :    

        cmap = plt.get_cmap()
        # qs = ax.quiver(X.flatten(), Y.flatten(), Z.flatten(), U.flatten(), V.flatten(), W.flatten(), pivot='tail', length=arrow_scale ) 
        # qs.set_color(s_m.to_rgba(amplitudes.flatten())) #set_array(np.transpose(amplitudes.flatten()))

        Us, Vs, Ws = disp_projected * arrow_scale / disp_projected.max()
        ax.scatter((X+Us).flatten(), (Y+Vs).flatten(), (Z+Ws).flatten(), c=s_m.to_rgba(amplitudes.flatten()), marker='.', linewidths=0)

        for i, x_val in np.ndenumerate(X):
            ax.plot([X[i], X[i]+Us[i]], [Y[i], Y[i]+Vs[i]], [Z[i], Z[i]+Ws[i]],'-', color=s_m.to_rgba(amplitudes[i]) , linewidth=1)
        

    ## For a color-coded surface representing amplitudes
    if style in ('*', 's', 'surf', 'surface'):

        ax.plot_surface(X*amplitudes_surf, Y*amplitudes_surf, Z*amplitudes_surf,linewidth=0.5, rstride=1, cstride=1, facecolors=s_m.to_rgba(amplitudes))    

    

    plt.show() 
    return ax


def energy_seismicsourcemodel(G, XYZ) :    
    """
    Evaluate statistical properties of the given seismic wave radiation 
    pattern.
    ______________________________________________________________________
    :type G : 3D array, list | np.array
    :param G : The vector of cartessian coordinates of the radiation 
        pattern.

    :type XYZ : 3D array, list | np.array
    :param XYZ : The cartessian coordinates of the origin points of the
        radiation pattern.

    :rtype : list
    :return : The statistical properties of amplitudes [rms, euclidian norm, average].

    .. todo::

        Compute more properties, add request handling.
    ______________________________________________________________________
    """

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
    """
    Set an instance of class Aki_Richards() that can be used for 
    earthquake modeling based on Aki and Richards (2002, eq. 4.29).
    ______________________________________________________________________
    :type mt : list
    :param mt : The focal mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - 
        the six independent components of the moment tensor).

    .. note::

        This object is composed of three methods : radpat, plot and 
        energy.

    .. seealso::

            :meth: `radpat`
            :meth: `plot`
            :meth: `energy`            
    ______________________________________________________________________

    """  

    def __init__(self, mt):
        self.mt = mt
        
    def radpat(self, wave='P', obs_cart='None', obs_sph='None'):
        """
        Returns the farfield radiation pattern (normalized displacement) based 
        on Aki and Richards (2002, eq. 4.29) and the cartesian coordinates of 
        the observation points.
        ______________________________________________________________________
        :type self : object: Aki_Richards
        :param self : This method use de self.mt attribute, i.e. the focal 
            mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - the six 
            independent components of the moment tensor)

        :type wave: String
        :param wave: type of wave to compute

        :type obs_cart : list | np.array 
        :param obs_cart : 3D vector array specifying the observations points
            in cartesian coordinates. 

        :type obs_sph : list | np.array 
        :param obs_sph : 3D vector array specifying the observations points
            in spherical coordinates (radians). The default is a unit sphere.

        :rtype : np.array 
        :return : 3D vector array with same shape than requested, that 
            contains the displacement vector for each observation point.

        .. rubric:: _`Supported wave`

            ``'P'``
                Bodywave P.

            ``'S' or 'S wave' or 'S-wave``
                Bodywave S.

            ``'Sv'``
                Projection of S on the meridian of the focal sphere.

            ``'Sh'``
                Projection of S on the parallels of the focal sphere. 

        .. note::

            This is based on MMesch, ObsPy, radpattern.py :
                https://github.com/MMesch/obspy/blob/radpattern/obspy/core/event/radpattern.py

        .. seealso::
    
            Aki, K., & Richards, P. G. (2002). Quantitative Seismology. (J. 
                Ellis, Ed.) (2nd ed.). University Science Books

        .. todo::

            Implement Sv and Sh wave (see below)
        ______________________________________________________________________
        """

        # Special cases ######################################################
        if wave in ('SH', 'Sv', 'S_v', 'Sv wave', 'Sv-wave'):

            ## Get S waves
            disp, obs_cart = self.radpat(wave='S', obs_cart=obs_cart, obs_sph=obs_sph)
            ## Project on Sv component
            disp = project_vectors(disp, vector_normal(obs_cart, 'v')) 

            return disp, obs_cart

        elif wave in ('SN', 'Sn', 'Snrh', 'Snrh wave', 'Snrh-wave'):

            ## Get S waves
            disp, obs_cart = self.radpat(wave='S', obs_cart=obs_cart, obs_sph=obs_sph)
            ## Project on Sh component
            disp = project_vectors(disp, vector_normal(obs_cart, 'nhr')) 

            return disp, obs_cart

        elif wave in ('SV', 'Sh', 'S_h', 'Sh wave', 'Sh-wave'):

            ## Get S waves
            disp, obs_cart = self.radpat(wave='S', obs_cart=obs_cart, obs_sph=obs_sph)
            ## Project on Sh component
            disp = project_vectors(disp, vector_normal(obs_cart, 'h')) 

            return disp, obs_cart
        ######################################################################

        # Get full mt
        Mpq = mt_full(self.mt)
        

        # Get observation points
        if obs_sph == 'None':
            obs_sph = sphere(r=1.,n=1000.)
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

                # loop through displacement component [n index]
                Mp = np.dot(Mpq, gamma)
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

    def plot(self, wave='P',style='*', comp=None, ax=None) :    
        """
        Plot the radiation pattern.
        ______________________________________________________________________
        :type self : object: Aki_Richards
        :param self : This method use de self.mt attribute, i.e. the focal 
            mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - the six 
            independent components of the moment tensor)

        :type wave: string
        :param wave: type of wave to compute

        :type style : string
        :param style : type of plot.

        :rtype : graphical object
        :return : The vector of cartessian coordinates.

        .. rubric:: _`Supported wave`

            ``'P'``
                Bodywave P.

            ``'S' or 'S wave' or 'S-wave``
                Bodywave S.

            ``'Sv'``
                Projection of S on the meridian of the focal sphere.

            ``'Sh'``
                Projection of S on the parallels of the focal sphere. 

        .. rubric:: _`Supported style`

            ``'None'``
                Plot the amplitudes as surface coding absolute amplitude as 
                radial distance from origin and amplitude polarity with color.
                 (uses :func:`numpy.abs`).

            ``'sign' or 'bb' or 'beachball``
                Plot a unit sphere, amplitude sign being coded with color.
                 (uses :func:`numpy.abs`).

            ``'frame'``
                Plot the amplitudes as a mesh, coding absolute amplitude as 
                radial distance. Amplitude polarity is not represented.

            ``'quiver'``
                Plot each vectors of the radiation pattern, sign and amplidude
                being coded with colors. 

        .. seealso::

                :func: `plot_seismicsourcemodel`
                :class: `Aki_Richards.radpat`
        ______________________________________________________________________

        """
        # Get radiation pattern
        G, XYZ = self.radpat(wave)
        if comp == None :
            if wave in ('Sv', 'Sv wave', 'Sv-wave'):
                comp='v'
            elif wave in ('Sn', 'Snrh','Snrh wave', 'Snrh-wave'):
                comp='n'
            elif wave in ('Sh', 'Sh wave', 'Sh-wave'):
                comp='h'
            elif wave in ('P', 'P wave', 'P-wave'):
                comp='r'

        plot_seismicsourcemodel(G, XYZ, style=style, mt=self.mt, comp=comp, ax=ax)

    def energy(self, wave='P') :   
        """
        Get statistical properties of radiation pattern.
        ______________________________________________________________________
        :type self : object: Aki_Richards
        :param self : This method use de self.mt attribute, i.e. the focal 
            mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - the six 
            independent components of the moment tensor)

        :type wave: string
        :param wave: type of wave to compute

        :rtype : list
        :return : The statistical properties of amplitudes [rms, euclidian norm, average].

        .. rubric:: _`Supported wave`

            ``'P'``
                Bodywave P.

            ``'S' or 'S wave' or 'S-wave``
                Bodywave S.

            ``'Sv'``
                Projection of S on the meridian of the focal sphere.

            ``'Sh'``
                Projection of S on the parallels of the focal sphere. 

        .. seealso::

                :func: `energy_seismicsourcemodel`
                :class: `Aki_Richards.radpat`
        ______________________________________________________________________

        """ 

        # Get radiation pattern and estimate energy
        G, XYZ = self.radpat(wave)  
        estimators = energy_seismicsourcemodel(G, XYZ)
        return estimators


class Vavryeuk(object):
    """
    Set an instance of class Vavryeuk() that can be used for 
    earthquake modeling based on Vavryèuk (2001).
    ______________________________________________________________________
    :type mt : list
    :param mt : The focal mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - 
        the six independent components of the moment tensor).

    :type poisson: variable
    :param poisson: Poisson coefficient (0.25 by default).

    .. note::

        This object is composed of three methods : radpat, plot and 
        energy.

    .. seealso::

            :meth: `radpat`
            :meth: `plot`
            :meth: `energy`            
    ______________________________________________________________________

    """  

    def __init__(self, mt, poisson=0.25):
        self.mt = mt
        self.poisson = poisson
        
    def radpat(self, wave='P', obs_cart='None', obs_sph='None'):
        """
        Returns the farfield radiation pattern (normalized displacement) based 
        on Vavryèuk (2001) and the cartesian coordinates of the observation 
        points.
        ______________________________________________________________________
        :type self : object: Vavryeuk
        :param self : This method use de self.poisson and self.mt attributes, 
            i.e. the focal mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - 
            the six independent components of the moment tensor).

        :type wave: String
        :param wave: type of wave to compute

        :type obs_cart : list | np.array 
        :param obs_cart : 3D vector array specifying the observations points
            in cartesian coordinates. 

        :type obs_sph : list | np.array 
        :param obs_sph : 3D vector array specifying the observations points
            in spherical coordinates (radians). The default is a unit sphere. 

        :rtype : np.array 
        :return : 3D vector array with same shape than requested, that 
            contains the displacement vector for each observation point.

        .. rubric:: _`Supported wave`

            ``'P'``
                Bodywave P.

            ``'S' or 'S wave' or 'S-wave``
                Bodywave S.

            ``'Sv'``
                Projection of S on the meridian of the focal sphere.

            ``'Sh'``
                Projection of S on the parallels of the focal sphere. 

        .. note::

            This is based on Kwiatek, G. (2013/09/15). Radiation pattern from 
            shear-tensile seismic source. Revision: 1.3 :
            http://www.npworks.com/matlabcentral/fileexchange/43524-radiation-pattern-from-shear-tensile-seismic-source

        .. seealso::            

            [1] Vavryèuk, V., 2001. Inversion for parameters of tensile 
             earthquakes.” J. Geophys. Res. 106 (B8): 16339–16355. 
             doi: 10.1029/2001JB000372.

            [2] Ou, G.-B., 2008, Seismological Studies for Tensile Faults. 
             Terrestrial, Atmospheric and Oceanic Sciences 19, 463.

            [3] Kwiatek, G. and Y. Ben-Zion (2013). Assessment of P and S wave 
             energy radiated from very small shear-tensile seismic events in 
             a deep South African mine. J. Geophys. Res. 118, 3630-3641, 
             doi: 10.1002/jgrb.50274

        .. rubric:: Example

            ex = NnK.source.SeismicSource([0,0,0,0,0,1])
            Svect, obs =  ex.Vavryeuk.radpat('S')                   

        .. plot::

            ex.Vavryeuk.plot()
        ______________________________________________________________________

        """
        # 1) Calculate moving rms and average, see moving()
        # 2) ...
        # 3) ... 
        
        poisson = self.poisson

        # Get angle from moment tensor
        [strike, dip, rake], [DC, CLVD, iso, deviatoric] = mt_angles(self.mt) 
        ## in radians
        [strike, dip, rake] = np.deg2rad([strike, dip, rake])
        ## convert DC ratio to angle
        MODE1 = np.arcsin((100-DC)/100.)        

        # Get observation points
        if obs_cart == 'None' :
            obs_cart = spherical_to_cartesian(sphere(r=1.,n=1000.))
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
            G = np.cos(TKO)*(np.cos(TKO)*(np.sin(MODE1)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(MODE1)*np.sin(rake)) - np.cos(AZM)*np.sin(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(MODE1)*np.sin(strike)) + np.sin(AZM)*np.sin(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(MODE1))) + np.sin(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(MODE1)) + np.cos(AZM)*np.sin(TKO)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)) + np.sin(AZM)*np.sin(TKO)*(np.cos(MODE1)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(AZM)*np.sin(TKO)*(np.cos(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(MODE1)*np.sin(strike)) - np.sin(AZM)*np.sin(TKO)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)) + np.cos(AZM)*np.sin(TKO)*(np.cos(MODE1)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
        
        elif wave in ('SH', 'Sh', 'Sh-wave', 'Sh wave', 'SH-wave', 'SH wave'):
            G = np.cos(TKO)*(np.cos(AZM)*(np.cos(MODE1)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(MODE1)) + np.sin(AZM)*(np.cos(MODE1)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(MODE1)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)) - np.cos(AZM)*(np.cos(MODE1)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(MODE1)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)))
            
        elif wave in ('SV', 'Sv', 'Sv-wave', 'Sv wave', 'SV-wave', 'SV wave'):
            G = np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)) - np.sin(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(MODE1)) + np.cos(TKO)*np.sin(AZM)*(np.cos(MODE1)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(MODE1)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(MODE1)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(MODE1)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(MODE1)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(MODE1))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(MODE1)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)) - np.cos(AZM)*np.cos(TKO)*(np.cos(MODE1)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))
        
        ## Re-using the same programme to get other things ... 
        elif wave in ('S', 'S-wave', 'S wave'):

            # for such definition this the less ugly
            Gsh = np.cos(TKO)*(np.cos(AZM)*(np.cos(MODE1)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(MODE1)) + np.sin(AZM)*(np.cos(MODE1)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(MODE1)*np.sin(strike))) - np.sin(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)) - np.cos(AZM)*(np.cos(MODE1)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) + np.cos(AZM)*np.sin(TKO)*(np.sin(AZM)*(np.cos(MODE1)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)) + np.cos(AZM)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)))
            Gsv = np.sin(AZM)*np.sin(TKO)*(np.cos(AZM)*np.cos(TKO)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)) - np.sin(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(MODE1)) + np.cos(TKO)*np.sin(AZM)*(np.cos(MODE1)*(np.sin(2*strike)*np.cos(rake)*np.sin(dip) - np.sin(2*dip)*np.cos(strike)**2*np.sin(rake)) - np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.cos(strike)**2*np.sin(dip)**2))) - np.cos(TKO)*(np.sin(TKO)*(np.sin(MODE1)*(2*np.cos(dip)**2 - (2*poisson)/(2*poisson - 1)) + np.sin(2*dip)*np.cos(MODE1)*np.sin(rake)) + np.cos(AZM)*np.cos(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(MODE1)*np.sin(strike)) - np.cos(TKO)*np.sin(AZM)*(np.cos(MODE1)*(np.cos(2*dip)*np.cos(strike)*np.sin(rake) - np.cos(dip)*np.cos(rake)*np.sin(strike)) - np.sin(2*dip)*np.cos(strike)*np.sin(MODE1))) + np.cos(AZM)*np.sin(TKO)*(np.sin(TKO)*(np.cos(MODE1)*(np.cos(2*dip)*np.sin(rake)*np.sin(strike) + np.cos(dip)*np.cos(rake)*np.cos(strike)) - np.sin(2*dip)*np.sin(MODE1)*np.sin(strike)) + np.cos(TKO)*np.sin(AZM)*(np.cos(MODE1)*(np.cos(2*strike)*np.cos(rake)*np.sin(dip) + (np.sin(2*dip)*np.sin(2*strike)*np.sin(rake))/2) - np.sin(2*strike)*np.sin(dip)**2*np.sin(MODE1)) - np.cos(AZM)*np.cos(TKO)*(np.cos(MODE1)*(np.sin(2*dip)*np.sin(rake)*np.sin(strike)**2 + np.sin(2*strike)*np.cos(rake)*np.sin(dip)) + np.sin(MODE1)*((2*poisson)/(2*poisson - 1) - 2*np.sin(dip)**2*np.sin(strike)**2)))

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


    def plot(self, wave='P',style='*', ax=None) :     
        """
        Plot the radiation pattern.
        ______________________________________________________________________
        :type self : object: Vavryeuk
        :param self : This method use de self.mt attribute, i.e. the focal 
            mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - the six 
            independent components of the moment tensor)

        :type wave: string
        :param wave: type of wave to compute

        :type style : string
        :param style : type of plot.

        :rtype : graphical object
        :return : The vector of cartessian coordinates.

        .. rubric:: _`Supported wave`

            ``'P'``
                Bodywave P.

            ``'S' or 'S wave' or 'S-wave``
                Bodywave S.

            ``'Sv'``
                Projection of S on the meridian of the focal sphere.

            ``'Sh'``
                Projection of S on the parallels of the focal sphere. 

        .. rubric:: _`Supported style`

            ``'None'``
                Plot the amplitudes as surface coding absolute amplitude as 
                radial distance from origin and amplitude polarity with color.
                 (uses :func:`numpy.abs`).

            ``'sign' or 'bb' or 'beachball``
                Plot a unit sphere, amplitude sign being coded with color.
                 (uses :func:`numpy.abs`).

            ``'frame'``
                Plot the amplitudes as a mesh, coding absolute amplitude as 
                radial distance. Amplitude polarity is not represented.

            ``'quiver'``
                Plot each vectors of the radiation pattern, sign and amplidude
                being coded with colors. 

        .. seealso::

                :func: `plot_seismicsourcemodel`
                :class: `Vavryeuk.radpat`
        ______________________________________________________________________

        """

        # Get radiation pattern and plot
        G, XYZ = self.radpat(wave)
        plot_seismicsourcemodel(G, XYZ, style=style, mt=self.mt, comp='r',ax=ax)

    def energy(self, wave='P') :   
        """
        Get statistical properties of radiation pattern.
        ______________________________________________________________________
        :type self : object: Vavryeuk
        :param self : This method use de self.mt attribute, i.e. the focal 
            mechanism NM x 6 (Mxx, Myy, Mzz, Mxy, Mxz, Myz - the six 
            independent components of the moment tensor)

        :type wave: string
        :param wave: type of wave to compute

        :rtype : list
        :return : The statistical properties of amplitudes [rms, euclidian norm, average].

        .. rubric:: _`Supported wave`

            ``'P'``
                Bodywave P.

            ``'S' or 'S wave' or 'S-wave``
                Bodywave S.

            ``'Sv'``
                Projection of S on the meridian of the focal sphere.

            ``'Sh'``
                Projection of S on the parallels of the focal sphere. 

        .. seealso::

                :func: `energy_seismicsourcemodel`
                :class: `Vavryeuk.radpat`
        ______________________________________________________________________

        """  
        # Get radiation pattern and estimate energy
        G, XYZ = self.radpat(wave)  
        estimators = energy_seismicsourcemodel(G, XYZ)
        return estimators



class SeismicSource(object):
    """
    Set an instance of class SeismicSource() that can be used for 
    earthquake modeling.
    ______________________________________________________________________
    :type mt : list
    :param mt : Definition of the focal mechanism supported 
        obspy.imaging.scripts.mopad.MomentTensor().

    :type poisson: variable
    :param poisson: Poisson coefficient (0.25 by default).

    .. example::

        ex = NnK.source.SeismicSource([1,2,3,4,5,6])
        ex.Aki_Richards.plot('P')

    .. note::

        This object is composed of three classes : MomentTensor, 
        Aki_Richards and Vavryeuk.

    .. seealso::

            :class: `obspy.imaging.scripts.mopad.MomentTensor`
            :class: `Aki_Richards`
            :class: `Vavryeuk`            
    ______________________________________________________________________

    """  
    # General attribut definition
    notes = 'Ceci est à moi'

    def __init__(self, mt, poisson=0.25):

        self.MomentTensor = MomentTensor(mt)
        self.Aki_Richards = Aki_Richards(np.asarray(self.MomentTensor.get_M('XYZ')))  
        self.Vavryeuk = Vavryeuk(np.asarray(self.MomentTensor.get_M('XYZ')),poisson = poisson)  


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





        
