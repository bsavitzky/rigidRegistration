"""
Save functions for stackregistration.py
"""

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
from os.path import splitext
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
import tifffile
import json

def save(imstack, fout, crop=True):
    """
    Saves imstack.average_image to fout.
    Saves as a 32-bit tif using tifffile package.

    Inputs:
        fout    str     path to output filename.
                        If fout does not end in .tif, it is appended
    """
    if splitext(fout)[1]=='.tif':
        filepath=fout
    else:
        filepath=fout+'.tif'
    metadata=dict()
    metadata['X_ij']=imstack.X_ij.tolist()
    metadata['Y_ij']=imstack.Y_ij.tolist()
    metadata['Rij_mask']=imstack.Rij_mask.tolist()
    metadata['shifts_x']=imstack.shifts_x.tolist()
    metadata['shifts_y']=imstack.shifts_y.tolist()
    metadata['xmin']=imstack.xmin
    metadata['xmax']=imstack.xmax
    metadata['ymin']=imstack.ymin
    metadata['ymax']=imstack.ymax
    metadata['method_of_cross_correlation']=imstack.correlation_type
    metadata['method_of_finding_maxima']=imstack.find_maxima_method
    if imstack.find_maxima_method=="gf":
        metadata['gaussian_num_of_peaks']=imstack.num_peaks
        metadata['gaussian_sigma_guess']=imstack.sigma_guess
        metadata['gaussian_window_radius']=imstack.window_radius
    metadata['fourier_mask_parameters']=imstack.mask_params
    metadata = json.dumps([metadata])
    if crop:
        tifffile.imsave(filepath,imstack.cropped_image.astype('float32'),description=metadata)
    else:
        tifffile.imsave(filepath,imstack.average_image.astype('float32'),description=metadata)
    return

def save_report(imstack, fout):

    # Set up report
    if splitext(fout)[1]=='.pdf':
        filepath=fout
    else:
        filepath=fout+'.pdf'
    report = PdfPages(fout)

    # Page 1: Image and FFT
    fig1,(ax11, ax12) = plt.subplots(1,2)
    ax11.matshow(imstack.average_image,cmap='gray')
    ax12.matshow(np.log(np.abs(np.fft.fftshift(np.fft.fft2(imstack.average_image))))[int(imstack.ny/4):int(3*imstack.ny/4),int(imstack.nx/4):int(3*imstack.nx/4)],cmap='gray',vmin=np.mean(np.log(np.abs(np.fft.fft2(imstack.average_image))).ravel()))
    ax11.axis('off')
    ax12.axis('off')
    fig1.tight_layout()
    fig1.suptitle("Average image")
    report.savefig()
    plt.close()

    # Page 2: Rij maps and mask

    # Make mask colormap
    cmap_mask=plt.cm.binary_r
    cmap_mask._init()
    alphas=np.linspace(1, 0, cmap_mask.N+3)
    cmap_mask._lut[:,-1] = alphas

    # Make figure
    fig2,((ax21,ax22),(ax23,ax24)) = plt.subplots(2,2)

    ax21.matshow(imstack.X_ij,cmap=r'RdBu')
    ax23.matshow(imstack.X_ij,cmap=r'RdBu')
    ax21.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax23.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    if np.sum(imstack.Rij_mask==False)!=0:
        ax23.matshow(imstack.Rij_mask,cmap=cmap_mask)

    ax22.matshow(imstack.Y_ij,cmap=r'RdBu')
    ax24.matshow(imstack.Y_ij,cmap=r'RdBu')
    ax22.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax24.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    if np.sum(imstack.Rij_mask==False)!=0:
        ax24.matshow(imstack.Rij_mask,cmap=cmap_mask)

    ax21.axis('off')
    ax22.axis('off')
    ax23.axis('off')
    ax24.axis('off')
    fig2.tight_layout()
    fig2.suptitle("Shift matrices")
    report.savefig()
    plt.close()

    # Page 3: Rij maps and mask

    # Make figure
    fig3,(ax31,ax32) = plt.subplots(1,2)

    ax31.matshow(imstack.X_ij_c,cmap=r'RdBu')
    ax31.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))

    ax32.matshow(imstack.Y_ij_c,cmap=r'RdBu')
    ax32.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))

    ax31.axis('off')
    ax32.axis('off')
    fig3.tight_layout()
    fig3.suptitle("Corrected shift matrices")
    report.savefig()
    plt.close()

    report.close()
    return


