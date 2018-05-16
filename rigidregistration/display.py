"""
Display functions for stackregistration.py
"""

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def show(imstack,crop=True,returnfig=False):
    """
    Show average image and its FFT.
    """
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(5,2.7),dpi=100)
    if crop:
        ax1.matshow(imstack.cropped_image,cmap='gray')
    else:
        ax1.matshow(imstack.average_image,cmap='gray')
    ax2.matshow(np.log(np.abs(np.fft.fftshift(np.fft.fft2(imstack.average_image))))[:imstack.nx,:imstack.ny],cmap='gray',vmin=np.mean(np.log(np.abs(np.fft.fft2(imstack.average_image))).ravel()))
    ax1.axis('off')
    ax2.axis('off')
    ax1.set_title("Averaged Image",y=1)
    ax2.set_title("Fourier Transform",y=1)
    ax1.grid(False)
    ax2.grid(False)
    fig.tight_layout()
    plt.subplots_adjust(left=0.03,right=0.97,bottom=0.03,top=0.89,wspace=0.03,hspace=0.01)
    if returnfig:
        return fig
    else:
        return

def show_Rij(imstack,Xmax=False,Ymax=False, mask=True,normalization=True,returnfig=False):
    """
    Display Rij matrix.

    Inputs:
        Xmax    float   Scales Xij colormap between -Xmax and +Xmax
        Ymax    float   Scales Yij colormap between -Ymax and +Ymax
        mask    bool    If true, overlays mask of bad data points.
    """

    X_ij_copy=np.copy(imstack.X_ij)
    Y_ij_copy=np.copy(imstack.Y_ij)
    ismask = True
    for i in range(imstack.nz):
        for j in range(imstack.nz):
            if imstack.Rij_mask[i,j]==False:
                ismask = False
                if normalization:
                    X_ij_copy[i,j]=0
                    Y_ij_copy[i,j]=0
    if mask and not (ismask and imstack.nz_min==0 and imstack.nz_max==imstack.nz ):

        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(5,2.7),dpi=100)
        if Xmax:
            ax1.matshow(X_ij_copy,cmap=r'RdBu',vmin=-Xmax,vmax=Xmax)
        else:
            ax1.matshow(X_ij_copy,cmap=r'RdBu')
        if Ymax:
            ax2.matshow(Y_ij_copy,cmap=r'RdBu',vmin=-Ymax,vmax=Ymax)
        else:
            ax2.matshow(Y_ij_copy,cmap=r'RdBu')

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

    ax1.grid(False)
    ax2.grid(False)
    ax1.set_title("Shift Matrix (X)",y=1.09)
    ax2.set_title("Shift Matrix (Y)",y=1.09)
    ax1.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax2.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax1.xaxis.set_ticks(np.arange(0, imstack.nz, 5))
    ax2.xaxis.set_ticks(np.arange(0, imstack.nz, 5))
    ax1.yaxis.set_ticks(np.arange(0, imstack.nz, 5))
    ax2.yaxis.set_ticks(np.arange(0, imstack.nz, 5))
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.03)
    if returnfig:
        return fig
    else:
        return


def show_Rij_c(imstack,Xmax=False,Ymax=False, mask=True):
    """
    Display corrected Rij matrix.

    Inputs:
        Xmax    float   Scales Xij colormap between -Xmax and +Xmax
        Ymax    float   Scales Yij colormap between -Ymax and +Ymax
        mask    bool    If true, overlays mask of bad data points.
    """
    fig,(ax1,ax2)=plt.subplots(1,2)
    if Xmax:
        ax1.matshow(imstack.X_ij_c,cmap=r'RdBu',vmin=-Xmax,vmax=Xmax)
    else:
        ax1.matshow(imstack.X_ij_c,cmap=r'RdBu')
    if Ymax:
        ax2.matshow(imstack.Y_ij_c,cmap=r'RdBu',vmin=-Ymax,vmax=Ymax)
    else:
        ax2.matshow(imstack.Y_ij_c,cmap=r'RdBu')
    if mask and np.sum(imstack.Rij_mask_c==False)!=0:
        # Make transparent colormap
        cmap_mask=plt.cm.binary_r
        cmap_mask._init()
        alphas=np.linspace(1, 0, cmap_mask.N+3)
        cmap_mask._lut[:,-1] = alphas
        # Overlay mask
        ax1.matshow(imstack.Rij_mask_c,cmap=cmap_mask)
        ax2.matshow(imstack.Rij_mask_c,cmap=cmap_mask)
    ax1.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax2.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))
    ax1.grid(False)
    ax2.grid(False)
    ax1.set_title("X shifts")
    ax2.set_title("Y shifts")
    plt.tight_layout()
    plt.show()
    return

def show_Fourier_mask(imstack,i=0,j=1):
    """
    Shows the mask used on cross correlations in Fourier space, overlaid on the Fourier
    transform of one image, and the cross correlation generated with this mask

    Inputs:
        i,j      ints     Image indices.  FFT displayed is of image i, cross correlation
                          displayed is between images i and j.
    """
    fig,(ax1,ax2,ax3)=plt.subplots(1,3)
    fig.suptitle(", ".join(["{} = {}".format(key,imstack.mask_params[key]) for key in list(imstack.mask_params)]))
    ax1.matshow(np.log(np.abs(np.fft.fftshift(imstack.fftstack[:,:,i]))),
                cmap='gray',vmin=np.average(np.log(np.abs(imstack.fftstack[:,:,i]))))
    ax1.matshow(np.fft.fftshift(imstack.mask_fourierspace),cmap='hot',alpha=0.4)
    ax2.matshow(np.log(np.abs(np.fft.fftshift(imstack.fftstack[:,:,i]*np.where(imstack.mask_fourierspace>0.0001,imstack.mask_fourierspace,0.0001)))), cmap='gray',
                vmin=1*np.average(np.log(np.abs(imstack.fftstack[:,:,i]))), vmax=1.8*np.average(np.log(np.abs(imstack.fftstack[:,:,i]))))
    ax3.matshow(np.abs(np.fft.fftshift(np.fft.ifft2(imstack.mask_fourierspace*imstack.fftstack[:,:,i]*imstack.fftstack[:,:,j]))),cmap='viridis')
    ax1.axis('off')
    ax2.axis('off')
    ax3.axis('off')
    ax1.grid(False)
    ax2.grid(False)
    ax3.grid(False)
    ax1.set_title("FFT with mask overlay")
    ax2.set_title("Masked FFT")
    ax3.set_title("Cross correlation")
    plt.show()
    return


def show_report(imstack):

    # Fig 1: Image and FFT
    fig1,(ax11, ax12) = plt.subplots(1,2)
    ax11.matshow(imstack.average_image,cmap='gray')
    ax12.matshow(np.log(np.abs(np.fft.fftshift(np.fft.fft2(imstack.average_image))))[int(imstack.ny/4):int(3*imstack.ny/4),int(imstack.nx/4):int(3*imstack.nx/4)],cmap='gray',vmin=np.mean(np.log(np.abs(np.fft.fft2(imstack.average_image))).ravel()))
    ax11.axis('off')
    ax12.axis('off')
    ax11.grid(False)
    ax12.grid(False)
    fig1.tight_layout()
    fig1.suptitle("Average image")
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
    ax21.grid(False)
    ax22.grid(False)
    ax23.grid(False)
    ax24.grid(False)
    fig2.tight_layout()
    fig2.suptitle("Shift matrices")
    plt.show()

    # Page 3: corrected Rij maps and mask

    # Make figure
    fig3,(ax31,ax32) = plt.subplots(1,2)

    ax31.matshow(imstack.X_ij_c,cmap=r'RdBu')
    ax31.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))

    ax32.matshow(imstack.Y_ij_c,cmap=r'RdBu')
    ax32.add_patch(Rectangle((imstack.nz_min-0.5, imstack.nz_min-0.5),imstack.nz_max-imstack.nz_min,imstack.nz_max-imstack.nz_min,facecolor='none',edgecolor='k',linewidth=3))

    ax31.axis('off')
    ax32.axis('off')
    ax31.grid(False)
    ax32.grid(False)
    fig3.tight_layout()
    fig3.suptitle("Corrected shift matrices")
    plt.show()

    return


