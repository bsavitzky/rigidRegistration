# -*- coding: utf-8 -*-
"""
Image registration optimized for noisy STEM datasets.

Designed for interactive use in the iPython notebook / Jupyter.
See, e.g., stackregistration_sample_notebook.ipynb.
"""

# Import global libraries
from __future__ import print_function, division, absolute_import
import numpy as np
from math import floor, ceil

# Import local libraries
from . import display
from . import save
from . import FFTW
from .utils import generateShiftedImage, gauss2d, fit_gaussian, on_edge, get_cutout, fit_peaks, getpaths, allpaths

class imstack(object):
    """
    Object for functions on stacks of images.
    Optimized to register and average stacks of noisy STEM images.
    """

    def __init__(self, image_stack):
        """
        Initializes imstack object from a 3D numpy array.

        Inputs:
            image_stack     ndarray of floats, shape (nx, ny, nz)
                            nz steps between images
        """
        self.imstack = image_stack
        self.nx, self.ny, self.nz = np.shape(image_stack)
        self.nz_min, self.nz_max = 0, self.nz

        # Define real and reciprocal space meshgrids
        rx,ry = np.meshgrid(np.arange(self.nx),np.arange(self.ny))
        self.rx,self.ry = rx.T,ry.T
        nx,ny = float(self.nx),float(self.ny)
        kx,ky = np.meshgrid(np.arange(-int(nx/2),int(nx/2),1),np.arange(-int(ny/2),int(ny/2),1))
        self.kx,self.ky = kx.T,ky.T
        self.kr = np.sqrt(self.kx**2+self.ky**2)

        # Define mask used to throw out bad datapoints
        self.nz_mask = np.ones((self.nz,self.nz),dtype=bool)
        self.bad_image_mask = np.ones((self.nz,self.nz),dtype=bool)
        self.outlier_mask = np.ones((self.nz,self.nz),dtype=bool)
        self.Rij_mask = np.ones((self.nz,self.nz),dtype=bool)

        # Include fftw library, if available
        self.fftw = FFTW.WrapFFTW((self.nx,self.ny))

        # Set initial params for gaussian fitting and CoM detection
        self.setGaussianFitParams()
        self.setCoMParams()

        return

    def getFFTs(self):
        """
        Gets all Fourier transforms of imstack.
        """
        mask_realspace = (np.sin(np.pi*self.rx/self.nx)*np.sin(np.pi*(self.ry/self.ny)))**2
        fftstack = np.zeros_like(self.imstack, dtype='complex')
        for i in range(self.nz):
            fftstack[:,:,i] = self.fftw.fft((self.imstack[:,:,i]-np.mean(self.imstack[:,:,i]))*mask_realspace)
        self.fftstack = fftstack
        return

    def makeFourierMask(self, mask='bandpass', n=4):
        """
        Defines the Fourier space mask to be used when finding cross-correlation.

        Inputs:
            n   int or float        Defines upper frequency allowed, in pixels - i.e. features
                                    smaller than ~n pixels are smoothed.  For masks with
                                    discrete cutoffs, the maximum frequency is k_max / n.
            mask    str             Either "bandpass", "lowpass", "hann", "hamming", "blackman",
                                    "gaussian", or "none".
                                    "lowpass" and "hann" are identical, returning a mask with a
                                    cos**2(k) functional form.
                                    "bandpass" uses a sin**2 functional form - i.e. decays to 0
                                    at k=0, but does not included a minimum cutoff frequency -
                                    bandpass masks with highpass cutoffs must be defined by the
                                    user.
                                    "gaussian" places the cutoff frequency at 3*sigma,
                                    containing 98.889% of the Fourier space intensity within
                                    k<k_max.
        """
        self.mask_params={'masktype':mask,'n':n}
        nx,ny = float(self.nx),float(self.ny)
        k_max = ny/n/2
        if mask=="bandpass":
            self.mask_fourierspace = np.fft.fftshift((self.kr<k_max)*(np.sin(2*n*np.pi*self.kr/ny)*np.sin(2*n*np.pi*self.kr/ny)))
        elif mask=="lowpass" or mask=="hann":
            self.mask_fourierspace = np.fft.fftshift((self.kr<k_max)*(np.cos(n*np.pi*self.kr/ny)**2))
        elif mask=="hamming":
            self.mask_fourierspace = np.fft.fftshift((self.kr<k_max)*((27./50.)+(23./50.)*np.cos((np.pi*self.kr)/k_max)))
        elif mask=="blackman":
            self.mask_fourierspace = np.fft.fftshift((self.kr<k_max)*((21./50.)+0.5*np.cos((np.pi*self.kr)/k_max)+(2./25.)*np.cos((2*np.pi*self.kr)/k_max)))
        elif mask=="gaussian":
            self.mask_fourierspace = np.fft.fftshift(np.exp(-(self.kr/(k_max/3.))**2))
        elif mask=="none":
            self.mask_fourierspace = np.ones_like(self.kr)
        else:
            print("Mask type must be 'bandpass', 'lowpass', or 'none'.")
            print("Alternatively, define a custom mask by setting the self.mask_fourierspace attribute manually.  The self.kr coordinates may be useful.")
            return
        return

    def makeFourierMask_eg(self,n1,n2,theta):
        """
        Defines an elliptical Gaussian Fourier space mask to be used when finding
        cross-correlation.

        Inputs:
            n1,n2 int or float  Defines upper frequencies allowed, in pixels - i.e. features
                                smaller than ~n pixels are smoothed - where n1 and n2 are
                                the mimumum feature sizes along the two primary axes.  Thus
                                The ellipticity is epsilon = n1/n2.
                                The cutoff frequency corresponding to the n's are set to 3*sigma
            theta   float       The angle of the n1 axis
        """
        self.mask_params={'masktype':'elliptical gaussian','n1':n1,'n2':n2,'theta':np.degrees(theta)}
        nx,ny = float(self.nx),float(self.ny)
        k1_max,k2_max = ny/n1/2.,ny/n2/2.
        sigma1,sigma2 = k1_max/3.,k2_max/3.

        a = (np.cos(theta)**2)/(2*sigma1**2) + (np.sin(theta)**2)/(2*sigma2**2)
        b = -(np.sin(2*theta))/(4*sigma1**2) + (np.sin(2*theta))/(4*sigma2**2)
        c = (np.sin(theta)**2)/(2*sigma1**2) + (np.cos(theta)**2)/(2*sigma2**2)

        self.mask_fourierspace = np.fft.fftshift(np.exp( -(a*self.kx**2 + 2*b*self.kx*self.ky + c*self.ky**2) ))
        return

    def makeUserDefinedFourierMask(self,mask):
        """
        Defines the Fourier space mask to be used when finding cross-correlation using a unique, user defined mask.
        This function expects a mask defined with the origin of k-space at the center of an array (i.e. at (fov/2,fov/2)).
        For masks defined with the origin of k-space at (0,0), the mask can be set directly with self.mask_fourierspace = mask.

        Inputs:
            mask    A mask, of shape (fov,fov).
        """
        self.mask_params={'masktype':"User defined"}
        self.mask_fourierspace = np.fft.fftshift(mask)
        return

    def findImageShifts(self, correlationType="cc", findMaxima="pixel", verbose=True):
        """
        Gets all image shifts.
        Proceeds as follows:
            (1) If fftstack has not been calculated, calculate it.  Note that if a Fourier
                mask has not been defined, the default bandpass filter is used.
            (2) Get correlation between image pair (i,j).  Use selected correlation type.
                Options listed below.
            (3) Find maximum in correlation.  Use selected maximum finding approach.  Options
                below.
            (4) Repeat 2-3 for all image pairs.  Store output shifts in Rij matrices.

        Inputs:
            correlationType     str     Correlation type.  Options:
                                        'cc' = cross correlation
                                        'mc' = mutual corraltion
                                        'pc' = phase correlation
            findMaxima          str     Method to find maximum value in correlations. Options:
                                        'pixel' = identify maximum pixel
                                        'gf' = fit a gaussian about n maximal pixels. Helps with
                                               unit cell jumps from sampling noise in atomic
                                               resolution data.  Change fitting parameters by
                                               calling self.setGaussianFitParams().
                                        'com' = center of mass. Finds center of mass twice, once
                                                over whole image, and again over small region
                                                about center from first iteration.
            verbose             bool    If True prints images being correlated.
        Outputs:
            X_ij, Y_ij    ndarrays of floats, shape (nz,nz), of calculated shifts in X and Y.
        """
        # Define shift matrices
        self.X_ij, self.Y_ij = np.zeros((self.nz,self.nz)), np.zeros((self.nz,self.nz))

        # If fftstack is not defined, get all FFTs
        if not hasattr(self,'fftstack'):
            self.getFFTs()

        # Define correlation function call
        if correlationType=="cc":
            getSingleCorrelation = self.getSingleCrossCorrelation
        elif correlationType=="mc":
            getSingleCorrelation = self.getSingleMutualCorrelation
        elif correlationType=="pc":
            getSingleCorrelation = self.getSinglePhaseCorrelation
        else:
            print("'correlationType' must be 'cc', 'mc', or 'pc'.")
            return
        self.correlation_type=correlationType
        self.find_maxima_method=findMaxima
        # Define maximum finder function call
        if findMaxima=="pixel":
            findMaxima = self.getSingleShift_pixel
        elif findMaxima=="gf":
            findMaxima = self.getSingleShift_gaussianFit
        elif findMaxima=="com":
            findMaxima = self.getSingleShift_com
        else:
            print("'findMaxima' must be 'pixel', 'gf', or 'com'.")
            return

        # Calculate all image shifts
        for i in range (0, self.nz-1):
            for j in range(i+1, self.nz):
                if verbose:
                    print("Correlating images {} and {}".format(i,j))
                cc = getSingleCorrelation(self.fftstack[:,:,i], self.fftstack[:,:,j])
                xshift, yshift = findMaxima(cc)
                if xshift<self.nx/2:
                        self.X_ij[i,j] = xshift
                else:
                        self.X_ij[i,j] = xshift-self.nx
                if yshift<self.ny/2:
                        self.Y_ij[i,j] = yshift
                else:
                        self.Y_ij[i,j] = yshift-self.ny

        # Fill in remaining skew-symmetric matrix elements
        for i in range (0, self.nz-1):
            for j in range(i+1, self.nz):
                self.X_ij[j,i] = -self.X_ij[i,j]
                self.Y_ij[j,i] = -self.Y_ij[i,j]

        return self.X_ij, self.Y_ij

    ########### Methods for correlating image pairs #############

    def getSingleCrossCorrelation(self, fft1, fft2):
        """
        Cross correlates two images from previously calculated ffts.
        Applies self.mask_fourierspace.  If undefined, masks using bandpass filter with n=4.
        (See self.makeFourierMask for more info.)
        """
        try:
            cross_correlation = np.abs(self.fftw.ifft(self.mask_fourierspace * fft2 * np.conj(fft1)))
        except AttributeError:
            self.makeFourierMask()
            cross_correlation = np.abs(self.fftw.ifft(self.mask_fourierspace * fft2 * np.conj(fft1)))
        return cross_correlation

    def getSingleMutualCorrelation(self, fft1, fft2):
        """
        Calculates mutual correlation function for two images from previously calculated ffts.
        Applies self.mask_fourierspace.  If undefined, masks using bandpass filter with n=4.
        (See self.makeFourierMask for more info.)
        """
        try:
            mutual_correlation = np.abs(self.fftw.ifft(self.mask_fourierspace * fft2 * np.conj(fft1) / np.sqrt(np.abs(fft2 * np.conj(fft1)))))
        except AttributeError:
            self.makeFourierMask()
            mutual_correlation = np.abs(self.fftw.ifft(self.mask_fourierspace * fft2 * np.conj(fft1) / np.sqrt(np.abs(fft2 * np.conj(fft1)))))
        return mutual_correlation

    def getSinglePhaseCorrelation(self, fft1, fft2):
        """
        Calculates phase correlation function for two images from previously calculated ffts.
        Applies self.mask_fourierspace.  If undefined, masks using bandpass filter with n=4.
        (See self.makeFourierMask for more info.)
        """
        try:
            phase_correlation = np.abs(self.fftw.ifft(self.mask_fourierspace * fft2 * np.conj(fft1) / (np.abs(fft2)**2) ))
        except AttributeError:
            self.makeFourierMask()
            phase_correlation = np.abs(self.fftw.ifft(self.mask_fourierspace * fft2 * np.conj(fft1) / (np.abs(fft2)**2) ))
        return phase_correlation

    ########### Methods for getting shifts from correlation maxima ########## 

    def getSingleShift_pixel(self, cc):
        """
        Calculates the shift between two images from their cross correlation by finding the
        maximum pixel.
        """
        xshift, yshift = np.unravel_index(np.argmax(cc),(self.nx,self.ny))
        return xshift, yshift

    def setGaussianFitParams(self,num_peaks=5,sigma_guess=2,window_radius=6):
        self.num_peaks=num_peaks
        self.sigma_guess=sigma_guess
        self.window_radius=window_radius
        return

    def getSingleShift_gaussianFit(self,cc):
        """
        Calculates the shift between two images from their cross correlation by fitting a
        gaussian to the cc maximum.
        Fits gaussians to the self.num_peaks maximum pixels and uses the point with the
        maximum *fit* rather than the maximum *pixel* intensity, to handle sampling issues.
        Alter fitting params with self.setGaussianFitParams(), or by manually setting the
        selt.num_peaks, self.sigma_guess, or self.window_radius attributes.

        Notes:
        (1) Gaussian fits make sense for atomic resolution data, but may not be appropriate
        elsewere, depending on your data.
        (2) Fitting multiple gaussians can handle sampling artifacts which can lead to
        identifying the incorrect maximum point in the cc.
        (3) Absent sampling problems, subpixel fitting for images stacks with ~10 or more
        images may not differ significantly from pixel fitting, as final shifts are calculated
        by averaging all the shifts to a given image.
        """
        all_shifts = self.get_n_cross_correlation_maxima(cc,self.num_peaks)

        data = np.fft.fftshift(cc)
        est_positions = all_shifts
        est_sigmas = np.ones_like(all_shifts)*self.sigma_guess
        est_params=[est_positions,est_sigmas]

        amplitudes, positions, sigmas, thetas, offsets, success_mask = fit_peaks(data,est_params,self.window_radius,print_mod=1, verbose=False)

        shift_x, shift_y = positions[np.argmax(offsets+amplitudes),:]
        return shift_x-np.shape(cc)[0]/2.0, shift_y-np.shape(cc)[1]/2.0

    def setCoMParams(self,num_iter=2,min_window_frac=3):
        self.num_iter=num_iter
        self.min_window_frac=min_window_frac
        return

    def getSingleShift_com(self,cc):
        """
        TODO: Document this function
        """
        ccs=np.fft.fftshift(cc)
        norm=np.sum(ccs)
        x_com, y_com = np.sum(ccs*self.rx)/norm, np.sum(ccs*self.ry)/norm

        # Iterate
        n_vals=np.linspace(self.min_window_frac,0,self.num_iter,endpoint=False)[::-1]
        for n in n_vals:
            r_com = np.sqrt((x_com-self.rx)**2+(y_com-self.ry)**2)
            weights = (r_com<self.ny/n/2)*(np.cos(n*np.pi*r_com/self.ny))**2
            ccs_weighted = ccs*weights
            norm = np.sum(ccs_weighted)
            x_com,y_com = np.sum(self.rx*ccs_weighted)/norm, np.sum(self.ry*ccs_weighted)/norm

        return x_com-self.nx/2.0, y_com-self.ny/2

    def get_n_cross_correlation_maxima(self,cc,n):
        """
        Gets the maximum n pixels in a cross correlation.
        """
        cc_shift=np.fft.fftshift(cc)
        shifts=np.zeros((n,2))
        for i in range(n):
            shifts[i,0],shifts[i,1]=np.unravel_index(np.argmax(cc_shift),np.shape(cc_shift))
            cc_shift[int(shifts[i,0]),int(shifts[i,1])]=0
        return shifts

    ########### Methods for masking Rij matrix #############

    def update_Rij_mask(self):
        """
        Rij_mask is comprised of:
            nz_mask: False outside the specified range of images between nz_min and nz_max
            bad_image_mask: False at specified bad images
            outlier_mask: False on datapoints determined to be outliers
        """
        self.Rij_mask = (self.nz_mask)*(self.bad_image_mask)*(self.outlier_mask)
        return

    def set_nz(self, nz_min, nz_max):
        """
        Sets range of images to include in averaging by setting self.nz_mask.

        Inputs:
            nz_min  int     first image in imstack to include
            nz_max  int     last image in imstack to include
        """
        self.nz_min, self.nz_max = nz_min, nz_max
        self.nz_mask = np.zeros((self.nz,self.nz),dtype=bool)
        self.nz_mask[nz_min:nz_max,nz_min:nz_max]=True
        self.update_Rij_mask()
        return

    def set_bad_images(self, bad_images):
        """
        Marks specified images as bad data, which won't be included in final average image.

        Inputs:
            bad_images      list of ints    indices of images to throw out
        """
        self.bad_images = list(bad_images)
        self.bad_image_mask = np.ones((self.nz,self.nz),dtype=bool)
        for image in bad_images:
            self.bad_image_mask[image,:]=False
            self.bad_image_mask[:,image]=False
        self.update_Rij_mask()
        return

    ############ Methods for outlier detection ###################

    def get_outliers(self, threshold, maxpaths=5):
        """
        Find outliers in Rij matrix, which will not be used in calculating final average image.
        Outliers are calculated by enforcing that correct transitivity holds between relative
        image shifts.

        Inputs:
            threshhold  float   Threshhold value controlling how much deviation from perfect
                                transitivity is consisdered acceptable.
            maxpaths    int     The number of transitivity relationships connecting two images
                                which are used to evaluate if a single image shift is correct
        """
        transitivity_scores=np.zeros_like(self.X_ij)
        for i in range(len(self.X_ij)-1):
            for j in range(i+1,len(self.X_ij)):
                paths = getpaths(i,j,maxpaths,self.nz)
                for p in paths:
                    pdx = np.array([self.X_ij[ip] for ip in p])
                    pdy = np.array([self.Y_ij[ip] for ip in p])
                    transitivity_scores[i,j] += np.sqrt((pdx.sum()-self.X_ij[j,i])**2+(pdy.sum()-self.Y_ij[j,i])**2)
        transitivity_scores /= maxpaths
        for i in range(len(self.X_ij)-1):
            for j in range(i+1,len(self.Y_ij)):
                transitivity_scores[j,i] = transitivity_scores[i,j]
        self.outlier_mask = transitivity_scores<threshold
        self.update_Rij_mask()
        return

    def get_outliers_NN(self, max_shift):
        """
        An alternative, simpler approach to outlier detection.
        Find outliers by looking at the difference between each shift matrix element and its
        nearest neighbor elements. If too many are larger than max_shift, the element is
        identified as an outlier.

        Inputs:
            max_shift   int or float    threshhold value for NN shifts
        """
        self.max_shift = max_shift

        x_mask = np.ones_like(self.X_ij,dtype=bool)
        y_mask = np.ones_like(self.Y_ij,dtype=bool)

        size=self.nz_max-self.nz_min
        for i in range(size):
            jump_max=3
            if i==0:
                imin=0
                jump_max=jump_max-1
            else:
                imin=i-1
            if i==size-1:
                imax=size-1
                jump_max=jump_max-1
            else:
                imax=i+2
            jump_max1=jump_max
            for j in range(i,size):
                jump_max=jump_max1
                if j==0:
                    jmin=0
                    jump_max=jump_max-1
                else:
                    jmin=j-1
                if j==size-1:
                    jmax=size-1
                    jump_max=jump_max-1
                else:
                    jmax=j+2
                x_jumps=0
                y_jumps=0
                for ia in range(imin,imax):
                    for ja in range(jmin,jmax):
                        if np.abs(self.X_ij[ia,ja]-self.X_ij[i,j])>float(max_shift):
                            x_jumps+=1
                        if np.abs(self.Y_ij[ia,ja]-self.Y_ij[i,j])>float(max_shift):
                            y_jumps+=1
                if x_jumps>jump_max:
                    x_mask[i,j]=False
                if y_jumps>jump_max:
                    y_mask[i,j]=False

        mask=x_mask*y_mask
        for i in range (0, self.nz_max-self.nz_min-1):
            for j in range(i+1, self.nz_max-self.nz_min):
                mask[j,i] = mask[i,j]
        self.outlier_mask = mask
        self.update_Rij_mask()
        return

    ############ Methods for reconstructing average image ############

    def make_corrected_Rij(self):
        maxpaths=5
        good_images=np.nonzero(np.all(self.Rij_mask==False,axis=1)==False)[0]
        temp_mask = np.copy(self.Rij_mask)
        self.X_ij_c,self.Y_ij_c = np.where(self.Rij_mask,self.X_ij,float('nan')),np.where(self.Rij_mask,self.Y_ij,float('nan'))
        count=1
        while np.all(temp_mask[good_images,:][:,good_images])==False:
            if count%1000==0:
                maxpaths *= 2
            for i in range(len(self.X_ij)-1):
                for j in range(i+1,len(self.Y_ij)):
                    if not temp_mask[i,j]:
                        n = 0.
                        x,y = 0.,0.
                        paths = getpaths(i,j,maxpaths,self.nz)
                        for p in paths:
                            if np.all([temp_mask[ip] for ip in p]):
                                x += np.array([self.X_ij_c[ip] for ip in p]).sum()
                                y += np.array([self.Y_ij_c[ip] for ip in p]).sum()
                                n += 1
                        if n:
                            self.X_ij_c[i,j],self.X_ij_c[j,i] = -x/n,x/n
                            self.Y_ij_c[i,j],self.Y_ij_c[j,i] = -y/n,y/n
            temp_mask = (np.isnan(self.X_ij_c)==False)*(np.isnan(self.Y_ij_c)==False)
            count += 1
        self.Rij_mask_c = temp_mask
        self.X_ij_c[np.isnan(self.X_ij_c)] = 0
        self.Y_ij_c[np.isnan(self.Y_ij_c)] = 0
        return

    def get_imshifts(self):
        """
        Get image shifts from Rij matrix by averaging rows.
        Ignores bad cross correlations, identified in self.Rij_mask; see self.update_Rij_mask
        for more details.
        Ignores images for which no good cross correlations are available.
        """
        try:
            self.shifts_x=np.sum(self.X_ij_c*self.Rij_mask_c,axis=1)/np.where(np.sum(self.Rij_mask_c,axis=1),np.sum(self.Rij_mask_c,axis=1),1)
            self.shifts_y=np.sum(self.Y_ij_c*self.Rij_mask_c,axis=1)/np.where(np.sum(self.Rij_mask_c,axis=1),np.sum(self.Rij_mask_c,axis=1),1)
        except AttributeError:
            self.make_corrected_Rij()
            self.shifts_x=np.sum(self.X_ij_c*self.Rij_mask_c,axis=1)/np.where(np.sum(self.Rij_mask_c,axis=1),np.sum(self.Rij_mask_c,axis=1),1)
            self.shifts_y=np.sum(self.Y_ij_c*self.Rij_mask_c,axis=1)/np.where(np.sum(self.Rij_mask_c,axis=1),np.sum(self.Rij_mask_c,axis=1),1)
        undetermined_shift_indices=[i for (i,x) in enumerate(np.isnan(self.shifts_x)*np.isnan(self.shifts_y)) if x]
        try:
            self.bad_images = list(set(self.bad_images)|set(undetermined_shift_indices))
        except AttributeError:
            self.bad_images = undetermined_shift_indices
        return

    def get_averaged_image(self, get_shifts=True, correct_Rij=True):
        """
        Calculates average image by shifting each image by previously calculated
        self.shifts_x/y.
        Shifted images are interpolated by shifting in Fourier space via the Fourier
        shift theorem.
        Stack of shifted images are stored in self.stack_registered.
        """
        if correct_Rij:
            self.make_corrected_Rij()
        if get_shifts:
            self.get_imshifts()

        good_image_indices = self.Rij_mask_c.sum(axis=1).nonzero()[0]
        self.stack_registered=np.zeros((self.nx,self.ny,len(good_image_indices)))
        for i in range(len(good_image_indices)):
            self.stack_registered[:,:,i]=generateShiftedImage(self.imstack[:,:,good_image_indices[i]],self.shifts_x[good_image_indices[i]],self.shifts_y[good_image_indices[i]])
        self.average_image = np.sum(self.stack_registered,axis=2)/float(len(good_image_indices))
        return

    def crop_image(self):
        """
        This function determines the min/max values in the final, averaged image which represent
        physically meaningful information, using the calculated shift values. The cropped image
        is then simply defined as the corresponding subset of the average image.
        """
        self.xmin,self.xmax = self.shifts_x.max(),self.shifts_x.min()
        self.ymin,self.ymax = self.shifts_y.max(),self.shifts_y.min()
        if self.xmax>0:
            self.xmax=self.nx
        if self.xmin<0:
            self.xmin=0
        if self.ymax>0:
            self.ymax=self.ny
        if self.ymin<0:
            self.ymin=0
        self.xmin=int(ceil(self.xmin))
        self.xmax=int(floor(self.xmax))
        self.ymin=int(ceil(self.ymin))
        self.ymax=int(floor(self.ymax))
        self.cropped_image=self.average_image[self.xmin:self.xmax,self.ymin:self.ymax]
        return

    ########################  Display methods #########################

    def show(self,crop=True, returnfig=False):
        """
        Show average image and its FFT.
        """
        if crop:
            self.crop_image()
        if returnfig:
            return display.show(self,crop=crop,returnfig=returnfig)
        else:
            display.show(self,crop=crop,returnfig=returnfig)
            return

    def show_Rij(self,Xmax=False,Ymax=False, mask=True,normalization=True, returnfig=False):
        """
        Display Rij matrix.

        Inputs:
            Xmax    float   Scales Xij colormap between -Xmax and +Xmax
            Ymax    float   Scales Yij colormap between -Ymax and +Ymax
            mask    bool    If true, overlays mask of bad data points.
        """
        if returnfig:
            return display.show_Rij(self,Xmax=Xmax,Ymax=Ymax,mask=mask,normalization=normalization,returnfig=returnfig)
        else:
            display.show_Rij(self,Xmax=Xmax,Ymax=Ymax,mask=mask,normalization=normalization,returnfig=returnfig)
            return

    def show_Rij_c(self,Xmax=False,Ymax=False, mask=True):
        """
        Display corrected Rij matrix.

        Inputs:
            Xmax    float   Scales Xij colormap between -Xmax and +Xmax
            Ymax    float   Scales Yij colormap between -Ymax and +Ymax
            mask    bool    If true, overlays mask of bad data points.
        """
        display.show_Rij_c(self,Xmax=Xmax,Ymax=Ymax,mask=True)
        return

    def show_Fourier_mask(self, i=0,j=1):
        """
        Shows the mask used on cross correlations in Fourier space, overlaid on the Fourier
        transform of one image.

        Inputs:
            image_index     int     FFT to display
        """
        return display.show_Fourier_mask(self,i=i,j=j)

    def show_report(self):
        """
        Displays a report showing the average image, its FFT, and all shifts with and without
        the mask used.
        """
        display.show_report(self)
        return

    ####################### Saving methods ######################

    def save(self,fout,crop=True):
        """
        Saves imstack.average_image to fout.
        Saves as a 32-bit tif using tifffile package.
        If 'crop' is True, crops output image.

        Inputs:
            fout    str     path to output filename.
                            If fout does not end in .tif, it is appended
            crop    bool    If True, returns copped image.
        """
        if crop:
            self.crop_image()
        save.save(self,fout=fout,crop=crop)
        return

    def save_report(self,fout):
        """
        Saves a report showing the average image, its FFT, and all shifts with and without
        the mask used.

        Inputs:
            fout    str     path to output filename.
                            If fout does not end in .pdf, it is appended
        """
        save.save_report(self,fout=fout)
        return

    #################### END IMSTACK OBJECT ####################




