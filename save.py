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

def save(imstack, fout):
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
    tifffile.imsave(filepath,imstack.average_image.astype('float32'))
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
    ax22.matshow(imstack.X_ij,cmap=r'RdBu')
    ax21.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax22.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    if np.sum(imstack.Rij_mask==False)!=0:
        ax22.matshow(imstack.Rij_mask,cmap=cmap_mask)

    ax23.matshow(imstack.Y_ij,cmap=r'RdBu')
    ax24.matshow(imstack.Y_ij,cmap=r'RdBu')
    ax23.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
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
    fig3,((ax31,ax32),(ax33,ax34)) = plt.subplots(2,2)

    ax31.matshow(imstack.X_ij_c,cmap=r'RdBu')
    ax32.matshow(imstack.X_ij_c,cmap=r'RdBu')
    ax31.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax32.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    if np.sum(imstack.Rij_mask_c==False)!=0:
        ax32.matshow(imstack.Rij_mask_c,cmap=cmap_mask)

    ax33.matshow(imstack.Y_ij_c,cmap=r'RdBu')
    ax34.matshow(imstack.Y_ij_c,cmap=r'RdBu')
    ax33.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax34.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    if np.sum(imstack.Rij_mask_c==False)!=0:
        ax34.matshow(imstack.Rij_mask_c,cmap=cmap_mask)

    ax31.axis('off')
    ax32.axis('off')
    ax33.axis('off')
    ax34.axis('off')
    fig3.tight_layout()
    fig3.suptitle("Corrected shift matrices")
    report.savefig()
    plt.close()

    report.close()
    return


