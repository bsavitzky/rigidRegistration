"""
Utility functions for stackregistration.py
"""

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import combinations, chain, islice


def generateShiftedImage(im, xshift, yshift):
    """
    Generates a shifted image, with supports for subpixel shifts.

    Inputs:
        im      ndarray, shape=(nx,ny)  Image to shift
        xshift  int or float            pixels to shift x
        yshift  int or float            pixels to shift y

    Outputs:
        shifted_im  ndarray, shape=(nx,ny)  Shifted image
    """

    # Define real space meshgrids
    nx, ny = np.shape(im)
    rxT,ryT = np.meshgrid(np.arange(nx),np.arange(ny))
    rx,ry = rxT.T,ryT.T
    nx,ny = float(nx),float(ny)

    w = -np.exp(-(2j*np.pi)*(xshift*rx/nx+yshift*ry/ny))
    shifted_im = np.abs(np.fft.ifft2(np.fft.ifftshift(np.fft.fftshift(np.fft.fft2(im))*w)))

    return shifted_im

# Define 2D Gaussian function
def gauss2d(xy_meshgrid, amplitude, x0, y0, sigma_x, sigma_y, theta, offset):
    # Returns result as a 1D array that can be passed to scipy.optimize.curve_fit
    (x,y) = xy_meshgrid
    x0, y0 = float(x0), float(y0)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + np.abs(amplitude)*np.exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) + c*((y-y0)**2) ))
    return g.ravel()


def fit_gaussian(x0, y0, sigma_x, sigma_y, theta, offset, data, plot=False, verbose=True):
    # Define mesh for input values and initial guess
    xy_meshgrid = np.meshgrid(range(np.shape(data)[0]),range(np.shape(data)[1]))
    initial_guess = (data[x0,y0], x0, y0, sigma_x, sigma_y, theta, offset)
    # Perform fit and pull out centers
    try:
        popt, pcov = curve_fit(gauss2d, xy_meshgrid, data.ravel(), p0=initial_guess)#, ftol=1.49012e-10, xtol=1.49012e-10)
    except RuntimeError:
        if verbose:
            print("Particle could not be fit to a 2D gaussian.  Returning guess parameters.")
        return np.array(initial_guess),None,False
    # Plotting for troubleshooting
    if plot:
        data_fitted = gauss2d(xs,ys, *popt)
        fig,ax=plt.subplots(1,1)
        ax.matshow(data,cmap='gray')
        ax.contour(xs,ys,data_fitted.reshape(np.shape(data)[0],np.shape(data)[1]),8, colors='w')
        plt.show()
    return popt, pcov, True


def on_edge(im,x0,y0,radius):
    x,y=np.shape(im)[0],np.shape(im)[1]
    if (x0-radius>=0 and y0-radius>=0 and x0+radius+1<x and y0+radius+1<y):
        return False
    else:
        return True

def get_cutout(im,x0,y0,radius):
    return im[int(x0-radius):int(x0+radius+1), int(y0-radius):int(y0+radius+1)]


def fit_peaks(data, est_params, window_radius, print_mod=100, verbose=True):
    """
    Inputs:
        data           ndarray, floats
        est_params     len 2 list of ndarrays, as follows:
        est_params[0]   -   (n,2) ndarray, positions (x,y)
        est_params[1]   -   (n,2) ndarray, sigmas (s_x,s_y)
        window_radius  int
    Outputs:
        fit_amplitudes    (n,) ndarray, floats
        fit_positions     (n,2) ndarray, floats, (x,y)
        fit_sigmas        (n,2) ndarray, floats, (s_x,s_y)
        fit_thetas        (n,) ndarray, floats
        fit_offsets       (n,) ndarray, floats
        fit_success_mask  (n,) ndarray, booleans
    """
    reference_fit_params = []
    reference_fit_success_mask = []
    for i in range(len(est_params[0])):
        #if i%print_mod==0:
        #    print "Fitting column {} of {}".format(i,len(est_params[0]))
        x0,y0,sigma_x,sigma_y = est_params[0][i,0],est_params[0][i,1],est_params[1][i,0],est_params[1][i,1]
        if not on_edge(data,x0,y0,window_radius):
            cutout = get_cutout(data,x0,y0,window_radius)
            popt, pcov, fit_success = fit_gaussian(window_radius,window_radius,sigma_x,sigma_y,0,0,cutout,plot=False,verbose=verbose)
            popt[1:3]=popt[1:3][::-1]+[x0,y0]-window_radius
            reference_fit_params.append(popt)
            reference_fit_success_mask.append(fit_success)
    reference_fit_params=np.array(reference_fit_params)
    reference_fit_success_mask=np.array(reference_fit_success_mask)

    fit_amplitudes=reference_fit_params[:,0]
    fit_positions=reference_fit_params[:,1:3]
    fit_sigmas=reference_fit_params[:,3:5]
    fit_thetas=reference_fit_params[:,5]
    fit_offsets=reference_fit_params[:,6]
    fit_success_mask=reference_fit_success_mask
    return fit_amplitudes, fit_positions, fit_sigmas, fit_thetas, fit_offsets, fit_success_mask

def makeslice(seq):
    """
    Produces a sequence of array elements from a sequence of integers
    i.e. [1,2,3] yields [(1,2),(2,3)]
    Inputs
        seq:    array_like of integers
    Returns
        slices: array_like of tuples
    """
    tups = []
    for i in range(len(seq)-1):
        tups += [(seq[i+1],seq[i])]
    return tups

def allpaths(i,j,maxpaths=200):
    """
    Finds all paths between integers i and j, returning as many as maxpaths.
    The number of paths grows as |j-i|!, so a cutoff is necessary for practicality.

    Inputs:
        i,j       ints		endpoints of path
        maxpaths  int		maximum number of paths to return
    Returns:
        index_paths			a list of length at most maxpaths, each element of which
							is a sequence which connects matrix elements i and j
    """
    if i>j:
        tmp = j
        j = i
        i = tmp
    if j-i < 2:
        return [[(j,i)]]
    n = range(i+1,j)
    combs = chain(*(combinations(n,l) for l in range(1,len(n)+1)[::1]))
    seq = [[i]+list(c)+[j] for c in islice(combs,maxpaths)]
    return map(makeslice, seq)

def getpaths(i,j,maxpaths,nz):
    """
    Finds a set of paths between integers i and j of length maxpaths.
    Selects paths in a sensible order, perferencing, in order (a) paths with elements between (i,j) (i.e. forward in time),
    and (b) shorter paths.

    Inputs:
        i,j       ints    endpoints of path
        maxpaths  int     maximum number of paths to return
    Returns:
        index_paths       a list of length at most maxpaths, each element of which is a sequence which connects
                          matrix elements i and j
    """
    if i>j:
        tmp = j
        j = i
        i = tmp
    if j-i<2:
        seq = [[i,j]]
    else:
        n = range(i+1,j)
        combs = chain(*(combinations(n,l) for l in range(1,len(n)+1)))
        seq = [[i]+list(c)+[j] for c in islice(combs,maxpaths)]

    a,b,count=0,0,0
    while(len(seq)<maxpaths):
        if count%2==0:
            if i-a>=0:
                a+=1
            else:
                b+=1
            count+=1
        else:
            if j+b<nz-1:
                b+=1
            else:
                a+=1
            count+=1
        i0,j0 = i-a,j+b

        n1,n2 = range(i-1,i0,-1),range(j0-1,i+1,-1)
        combs1 = chain(*(combinations(n1,l) for l in range(1,len(n1)+1)[::-1]))
        combs2 = chain(*(combinations(n2,l) for l in range(1,len(n2)+1)[::-1]))
        seq1 = [list(c) for c in islice(combs1,maxpaths)]
        seq2 = [list(c) for c in islice(combs2,maxpaths)]
        seq1.append([])
        seq2.append([])
        for seq11 in seq1:
            for seq22 in seq2:
                subseq=[]
                if i!=i0:
                    subseq.append(i0)
                if j!=j0:
                    subseq.append(j0)
                seq.append([i]+seq11+subseq+seq22+[j])

    return map(makeslice, seq[:maxpaths])










