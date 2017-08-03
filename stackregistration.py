# -*- coding: utf-8 -*-
"""
Image registration optimized for noisy STEM datasets.

Designed for interactive use in the iPython notebook / Jupyter.
See, e.g., stackregistration_sample_notebook.ipynb.
"""

# Import global libraries
from __future__ import print_function, division, absolute_import
import numpy as np
import FFTW

# Import local libraries
import display
import save
from utils import generateShiftedImage, gauss2d, fit_gaussian, on_edge, get_cutout, fit_peaks

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

        # Set initial params for gaussian fitting
        self.setGaussianFitParams()

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
            n   int or float        Defines upper frequency allowed.  For all masks,
                                    maximum frequency is k_max / n.  Thus features smaller
                                    than ~n pixels are smoothed.
                                    Heuristically n=4 works well, but will vary with data.
            mask    str             Options:
                                    "bandpass" - a sin^2 mask which goes to zero at k=0
                                    "lowpass" - a sin^2 mask with high frequency cutoff defined
                                                by n
                                    "blackman" - a blackman lowpass mask
                                    "gaussian" - a gaussian lowpass mask.  High freq. cutoff
                                        defined by n is NOT a hard cutoff, and is set to 3*sigma
                                    "none" - no mask
        """
        nx,ny = float(self.nx),float(self.ny)
        if mask=="bandpass":
            self.mask_fourierspace = np.fft.fftshift((self.kr<ny/n/2)*(np.sin(2*n*np.pi*self.kr/ny)*np.sin(2*n*np.pi*self.kr/ny)))
        elif mask=="lowpass":
            self.mask_fourierspace = np.fft.fftshift((self.kr<ny/n/2)*(np.cos(n*np.pi*self.kr/ny)*np.cos(n*np.pi*self.kr/ny)))
            #self.mask_fourierspace = np.fft.fftshift(np.exp(-(4*s.kr/s.ny)**2))
        elif mask=="blackman":
            a=self.ny/n/2
            self.mask_fourierspace = np.fft.fftshift((self.kr<a)*((21./50.)+0.5*np.cos((np.pi*self.kr)/a)+(2./25.)*np.cos((2*np.pi*self.kr)/a)))
        elif mask=="none":
            self.mask_fourierspace = np.ones_like(self.kr)
        elif mask=="gaussian":
            sigma=self.ny/n/6.
            self.mask_fourierspace = np.fft.fftshift(np.exp(-self.kr**2/(2*sigma**2)))
        else:
            print("Mask type must be 'bandpass', 'lowpass', 'blackman', 'gaussian', or 'none'.")
            print("Alternatively, define a custom mask by setting the self.mask_fourierspace attribute manually.  The self.kr coordinates may be useful.")
            return
        return

    def makeEllipticalGaussianMask(self, n, theta, epsilon):
        """
        Defines the Fourier space mask to be used when finding cross-correlations as a
        non-isotropic gaussian function.

        Inputs:
            n   int or float    Defines high frequency 'cutoff'.  Not strictly a cutoff -
                                the frequency corresponding to n pixels is set to 3*sigma along
                                the axis oriented along theta.
            theta   float       Defines orientation of axes, in degrees
            epsilon float       Anisotropy of mask (axis1length/axis2length)
        """
        sigma1 = self.ny/n/6.
        sigma2 = sigma1*epsilon
        thetarad = np.pi/2.0 - np.radians(theta)

        a = (np.cos(thetarad)**2)/(2*sigma1**2) + (np.sin(thetarad)**2)/(2*sigma2**2)
        b = -(np.sin(2*thetarad))/(4*sigma1**2) + (np.sin(2*thetarad))/(4*sigma2**2)
        c = (np.sin(thetarad)**2)/(2*sigma1**2) + (np.cos(thetarad)**2)/(2*sigma2**2)

        self.mask_fourierspace = np.fft.fftshift(np.exp( -(a*self.kx**2 + 2*b*self.kx*self.ky + c*self.ky**2) ))
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

        # Define maximum finder function call
        if findMaxima=="pixel":
            findMaxima = self.getSingleShift_pixel
        elif findMaxima=="gf":
            findMaxima = self.getSingleShift_gaussianFit
        else:
            print("'findMaxima' must be 'pixel', or 'gf'.")
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

    def get_outliers(self, method="NN", *args):
        """
        Find outliers in Rij matrix, which will not be used in calculating final average image.

        Inputs:
            method  str     Method to be used for outlier detection.
                            Currenly supported: 'NN'

        Currently supported outlier detection methods:
        NN - detects outliers by looking at the shifts of the nearest neighbor image pairs
        args:
            arg[0]: max_shift - outlier shift threshhold
        """
        if method=="NN":
            self.outlier_mask = self.get_outliers_NN(args[0])
        else:
            print("Outlier detection method must be 'NN'. Skipping outlier detection")
        self.update_Rij_mask()
        return


    ############ Methods for outlier detection ###################

    def get_outliers_NN(self, max_shift):
        """
        Find outliers by looking at shift value to all 8 NN pixels in Rij.  If too many
        are larger than max_shift, cc is identified as an outlier.

        Inputs:
            max_shift   int or float    threshhold value for NN shifts
        """
        self.max_shift = max_shift

        x_mask = np.ones_like(self.X_ij,dtype=bool)
        y_mask = np.ones_like(self.Y_ij,dtype=bool)

        for i in range(self.nz_max-self.nz_min):
            for j in range(i,self.nz_max-self.nz_min):
                # Corners - more than one jump
                if (i==0 and j==0):
                    num_jumps=np.sum(np.abs(self.X_ij[i:i+2,j:j+2]-self.X_ij[i,j])>max_shift)
                    if num_jumps>1:
                        x_mask[i,j]=False
                elif (i==0 and j==self.nz_max-self.nz_min-1):
                    num_jumps=np.sum(np.abs(self.X_ij[i:i+2,j-1:j+1]-self.X_ij[i,j])>max_shift)
                    if num_jumps>1:
                        x_mask[i,j]=False
                elif (i==self.nz_max-self.nz_min-1 and j==0):
                    num_jumps=np.sum(np.abs(self.X_ij[i-1:i+1,j:j+2]-self.X_ij[i,j])>max_shift)
                    if num_jumps>1:
                        x_mask[i,j]=False
                elif (i==self.nz_max-self.nz_min-1 and j==self.nz_max-self.nz_min-1):
                    num_jumps=np.sum(np.abs(self.X_ij[i-1:i+1,j-1:j+1]-self.X_ij[i,j])>max_shift)
                    if num_jumps>1:
                        x_mask[i,j]=False
                # Edges
                elif (i==0):
                    num_jumps=np.sum(np.abs(self.X_ij[i:i+2,j-1:j+2]-self.X_ij[i,j])>max_shift)
                    if num_jumps>2:
                        x_mask[i,j]=False
                elif (i==self.nz_max-self.nz_min-1):
                    num_jumps=np.sum(np.abs(self.X_ij[i-1:i+1,j-1:j+2]-self.X_ij[i,j])>max_shift)
                    if num_jumps>2:
                        x_mask[i,j]=False
                elif (j==0):
                    num_jumps=np.sum(np.abs(self.X_ij[i-1:i+2,j:j+2]-self.X_ij[i,j])>max_shift)
                    if num_jumps>2:
                        x_mask[i,j]=False
                elif (j==self.nz_max-self.nz_min-1):
                    num_jumps=np.sum(np.abs(self.X_ij[i-1:i+2,j-1:j+1]-self.X_ij[i,j])>max_shift)
                    if num_jumps>2:
                        x_mask[i,j]=False
                # Bulk
                else:
                    num_jumps=np.sum(np.abs(self.X_ij[i-1:i+2,j-1:j+2]-self.X_ij[i,j])>max_shift)
                    if num_jumps>3:
                        x_mask[i,j]=False

        for i in range(self.nz_max-self.nz_min):
            for j in range(i,self.nz_max-self.nz_min):
                # Corners - more than one jump
                if (i==self.nz_min and j==self.nz_min):
                    num_jumps=np.sum(np.abs(self.Y_ij[i:i+2,j:j+2]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>1:
                        y_mask[i,j]=False
                elif (i==self.nz_min and j==self.nz_max-1):
                    num_jumps=np.sum(np.abs(self.Y_ij[i:i+2,j-1:j+1]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>1:
                        y_mask[i,j]=False
                elif (i==self.nz_max-1 and j==self.nz_min):
                    num_jumps=np.sum(np.abs(self.Y_ij[i-1:i+1,j:j+2]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>1:
                        y_mask[i,j]=False
                elif (i==self.nz_max-1 and j==self.nz_max-1):
                    num_jumps=np.sum(np.abs(self.Y_ij[i-1:i+1,j-1:j+1]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>1:
                        y_mask[i,j]=False
                # Edges
                elif (i==self.nz_min):
                    num_jumps=np.sum(np.abs(self.Y_ij[i:i+2,j-1:j+2]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>2:
                        y_mask[i,j]=False
                elif (i==self.nz_max-1):
                    num_jumps=np.sum(np.abs(self.Y_ij[i-1:i+1,j-1:j+2]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>2:
                        y_mask[i,j]=False
                elif (j==self.nz_min):
                    num_jumps=np.sum(np.abs(self.Y_ij[i-1:i+2,j:j+2]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>2:
                        y_mask[i,j]=False
                elif (j==self.nz_max-1):
                    num_jumps=np.sum(np.abs(self.Y_ij[i-1:i+2,j-1:j+1]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>2:
                        y_mask[i,j]=False
                # Bulk
                else:
                    num_jumps=np.sum(np.abs(self.Y_ij[i-1:i+2,j-1:j+2]-self.Y_ij[i,j])>max_shift)
                    if num_jumps>3:
                        y_mask[i,j]=False

        mask=x_mask*y_mask
        for i in range (0, self.nz_max-self.nz_min-1):
            for j in range(i+1, self.nz_max-self.nz_min):
                mask[j,i] = mask[i,j]

        return mask


    ############ Methods for reconstructing average image ############

    def get_imshifts(self):
        """
        Get image shifts from previously calculated Rij matrix by averaging rows.
        Ignores bad cross correlations, indentified in self.Rij_mask; see self.update_Rij_mask
        for more info.
        Ignores images for which no good cross correlations are available.
        """
        self.shifts_x=(1.0/np.sum(self.Rij_mask,axis=1))*np.sum(self.X_ij*self.Rij_mask,axis=1)
        self.shifts_y=(1.0/np.sum(self.Rij_mask,axis=1))*np.sum(self.Y_ij*self.Rij_mask,axis=1)

        undetermined_shift_indices=[i for (i,x) in enumerate(np.isnan(self.shifts_x)*np.isnan(self.shifts_y)) if x]
        try:
            self.bad_images = list(set(self.bad_images)|set(undetermined_shift_indices))
        except AttributeError:
            self.bad_images = undetermined_shift_indices
        return


    def get_averaged_image(self, get_shifts=True):
        """
        Calculates average image by shifting each image by previously calculated
        self.shifts_x/y.
        Shifted images are interpolated by shifting in Fourier space via the Fourier
        shift theorem.
        Stack of shifted images are stored in self.stack_registered.
        """
        if get_shifts:
            self.get_imshifts()

        good_image_indices = [i for i in np.arange(self.nz_min,self.nz_max) if i not in self.bad_images]
        self.stack_registered=np.zeros((self.nx,self.ny,len(good_image_indices)))
        for i in range(len(good_image_indices)):
            self.stack_registered[:,:,i]=generateShiftedImage(self.imstack[:,:,good_image_indices[i]],self.shifts_x[good_image_indices[i]],self.shifts_y[good_image_indices[i]])
        self.average_image = np.sum(self.stack_registered,axis=2)/float(len(good_image_indices))
        return


########################  Display methods #########################

    def show(self):
        """
        Show average image and its FFT.
        """
        display.show(self)
        return

    def show_Rij(self,Xmax=False,Ymax=False, mask=True):
        """
        Display Rij matrix.

        Inputs:
            Xmax    float   Scales Xij colormap between -Xmax and +Xmax
            Ymax    float   Scales Yij colormap between -Ymax and +Ymax
            mask    bool    If true, overlays mask of bad data points.
        """
        display.show_Rij(self,Xmax=Xmax,Ymax=Ymax,mask=True)
        return

    def show_Fourier_mask(self, image_index=0):
        """
        Shows the mask used on cross correlations in Fourier space, overlaid on the Fourier
        transform of one image.

        Inputs:
            image_index     int     FFT to display
        """
        display.show_Fourier_mask(self,image_index=image_index)
        return

    def show_report(self):
        """
        Displays a report showing the average image, its FFT, and all shifts with and without
        the mask used.
        """
        display.show_report(self)
        return

####################### Saving methods ######################

    def save(self, fout):
        """
        Saves imstack.average_image to fout.
        Saves as a 32-bit tif using tifffile package.

        Inputs:
            fout    str     path to output filename.
                            If fout does not end in .tif, it is appended
        """
        save.save(self,fout=fout)
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




