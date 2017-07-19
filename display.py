"""
Display functions for stackregistration.py
"""

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def show(imstack,crop=False):
    """
    Show average image and its FFT.
    """
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(5,2.7),dpi=100)
    if crop:
        ax1.matshow(imstack.cropped_image,cmap='gray')
    else:
        ax1.matshow(imstack.average_image,cmap='gray')
    ax2.matshow(np.log(np.abs(np.fft.fftshift(np.fft.fft2(imstack.average_image))))[:int(imstack.nx),:imstack.ny],cmap='gray',vmin=np.mean(np.log(np.abs(np.fft.fft2(imstack.average_image))).ravel()))
    ax1.axis('off')
    ax2.axis('off')
    ax1.set_title("Averaged Image")
    ax2.set_title("Fourier Transform")
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.2,hspace=0.05)
    return fig

def show_Rij(imstack,Xmax=False,Ymax=False, mask=True,normalization=False):
    """
    Display Rij matrix.

    Inputs:
        Xmax    float   Scales Xij colormap between -Xmax and +Xmax
        Ymax    float   Scales Yij colormap between -Ymax and +Ymax
        mask    bool    If true, overlays mask of bad data points.
    """
    
    xthing=np.copy(imstack.X_ij)
    ything=np.copy(imstack.Y_ij)
    meh = True
    for i in range(imstack.nz):
        for j in range(imstack.nz):
            if imstack.Rij_mask[i,j]==False:
                meh = False
                if normalization:
                    xthing[i,j]=0
                    ything[i,j]=0
    if mask and not (meh and imstack.nz_min==0 and imstack.nz_max==imstack.nz ):
        
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(5,2.7),dpi=100)
        if Xmax:
            ax1.matshow(xthing,cmap=r'RdBu',vmin=-Xmax,vmax=Xmax)
        else:
            ax1.matshow(xthing,cmap=r'RdBu')
        if Ymax:
            ax2.matshow(ything,cmap=r'RdBu',vmin=-Ymax,vmax=Ymax)
        else:
            ax2.matshow(ything,cmap=r'RdBu')
        
        # Make transparent colormap
        cmap_mask=plt.cm.binary_r
        cmap_mask._init()
        alphas=np.linspace(1, 0, cmap_mask.N+3)
        cmap_mask._lut[:,-1] = alphas
        # Make mask with full size
        full_mask = np.zeros_like(imstack.X_ij,dtype=bool)
        full_mask=imstack.Rij_mask
        imstack.full_mask=full_mask
        # Overlay mask
        ax1.matshow(full_mask,cmap=cmap_mask)
        ax2.matshow(full_mask,cmap=cmap_mask)
        ax1.set_title("Masked Shift Mastrix, X",y=1.1)
        ax2.set_title("Masked Shift Matrix, Y",y=1.1)
    else:
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(5,2.7),dpi=100)
        if Xmax:
            ax1.matshow(imstack.X_ij,cmap=r'RdBu',vmin=-Xmax,vmax=Xmax)
        else:
            ax1.matshow(imstack.X_ij,cmap=r'RdBu')
        if Ymax:
            ax2.matshow(imstack.Y_ij,cmap=r'RdBu',vmin=-Ymax,vmax=Ymax)
        else:
            ax2.matshow(imstack.Y_ij,cmap=r'RdBu')
        ax1.set_title("Shift Matrix in X Direction",y=1.1)
        ax2.set_title("Shift Matrix in Y Direction",y=1.1)
    ax1.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax2.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax1.xaxis.set_ticks(np.arange(0, imstack.nz, 5))
    ax2.xaxis.set_ticks(np.arange(0, imstack.nz, 5))
    ax1.yaxis.set_ticks(np.arange(0, imstack.nz, 5))
    ax2.yaxis.set_ticks(np.arange(0, imstack.nz, 5))
    plt.tight_layout()
    return fig

def show_Fourier_mask(imstack, image_index=0):
    """
    Shows the mask used on cross correlations in Fourier space, overlaid on the Fourier
    transform of one image.

    Inputs:
        image_index     int     FFT to display
    """
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(5,2.7),sharex=True,sharey=True)
    ax1.matshow(np.log(np.abs(np.fft.fftshift(imstack.fftstack[:,:,image_index]))),
                cmap='gray',vmin=np.mean(np.log(np.abs(np.fft.fftshift(imstack.fftstack[:,:,image_index]))).ravel()))
    ax1.matshow(np.fft.fftshift(imstack.mask_fourierspace),cmap='hot',alpha=0.3)
    ax2.matshow(np.log(np.abs(np.fft.fftshift(imstack.fftstack[:,:,image_index]*np.where(imstack.mask_fourierspace,imstack.mask_fourierspace,0.001)))),
                cmap='gray',vmin=np.mean(np.log(np.abs(np.fft.fftshift(imstack.fftstack[:,:,image_index]))).ravel()))
    ax1.axis('off')
    ax2.axis('off')
    ax1.set_title("FFT with mask overlay")
    ax2.set_title("Masked FFT")
    return fig

def show_report(imstack):

    # Fig 1: Image and FFT
    fig1,(ax11, ax12) = plt.subplots(1,2)
    ax11.matshow(imstack.average_image,cmap='gray')
    ax12.matshow(np.log(np.abs(np.fft.fftshift(np.fft.fft2(imstack.average_image))))[int(imstack.ny/4):int(3*imstack.ny/4),int(imstack.nx/4):int(3*imstack.nx/4)],cmap='gray',vmin=np.mean(np.log(np.abs(np.fft.fft2(imstack.average_image))).ravel()))
    ax11.axis('off')
    ax12.axis('off')
    fig1.tight_layout()
    plt.show()

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
    plt.show()

    return


