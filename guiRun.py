#!/usr/bin/env python

"""
- Store defaults for calculations for users cause users are picky
- See if we can make things better compared to window size
- Make sure everything is idiot proof so you don't get errors
- Make it so registered images don't wrap around (wrap around, go for the personal glory...)
- Bind mouse wheel to scroll bar cause its so much nicer
- Add pictures to user guide
- Labels being off - Make all labels same size and stuff
- Take out com radio button
- Mask look at it in popup window with 1 cross correlation
- Make user guides more bullet points 
- If you zoom in average image, shows Fourier mask just for that part of image
- Write metadata of all the stuff we did 
"""

# Import global libraries
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from time import time
from PIL import Image
from os import listdir
import subprocess
from tifffile import imread
import math as math
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
else:
    import tkinter as tk
import ttk as ttk
from tkFileDialog import askopenfilename #ttk overwrites this one 
import tkFileDialog
import Tkconstants
import os.path
import tkMessageBox
import matplotlib.backends.backend_tkagg as tkagg
import warnings
warnings.filterwarnings("ignore") #Theres a warning that prints if you don't do this (Don't worry, the warning is not important)
import stackregistration as stackregister #This imports the local library
import serReader

#Initialize global variables (will be given values later)
s=None
rijMaskFalse=None
rijMaskTrue=None
averageImage=None
fourierMask=None
canvas=None
frame=None
outpop = None
nzpop = None
outCoeEntry = None
outZEntry = None
outShiftEntry = None
fourierpop = None

loaded=False
sliceUn=-1
sliceReg=-1

LARGE_FONT = ("Verdana",24)
MEDIUMLARGE = ("Verdana",20)
SMALLFONT=("Verdana",12)

outlierMethod="NN"
zscore=2.0
shiftMax=20
polycoeffmax=3
correlationType="cc"
findMaxima="pixel"
fouriern = 4
fourierMaskType="bandpass"
correlationType="cc"
findMaxima="pixel"
gaussiannumpeaks=3
sigmaguess=2
windowradius=6
numiter=3
minwindowfrac=3

# This makes the popup and the functions for making a custom mask 
class CustomMask():
    def __init__(self,master):
        global loaded
        if loaded:
            self.master = master
            self.pop = tk.Tk()
            self.i=0
            self.j=0
            self.img=s.show_Rij(mask=True,normalization=True)
            self.maskData=np.copy(s.Rij_mask)
            self.original=np.copy(s.Rij_mask)
            self.frame=tk.Frame(self.pop)
            self.run()
        else:
            tkMessageBox.showwarning("Load Data","Please load images in to analze first")
    def run(self):
        self.pop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
        self.pop.title('Custom Outlier Mask')
        self.pop.resizable(width=False, height=False)
        self.drawCanvas()
        self.frame.grid(column=0,row=0,columnspan=2)
        saveBtn=tk.Button(self.pop,text="Save",font=MEDIUMLARGE,command=self.saving)
        cnclBtn=tk.Button(self.pop,text="Cancel",font=MEDIUMLARGE,command=self.cancel)
        saveBtn.grid(column=0,row=1)
        cnclBtn.grid(column=1,row=1)
    def saving(self):
        s.Rij_mask=self.maskData
        self.pop.destroy()
        recalc()
    def cancel(self):
        self.pop.destroy()
    def OnMouseDown(self,event):
        if sys.platform == "win32":
            xmin1=35
            xmax1=314
            ymin=46
            ymax=324
            xmin2=352
            xmax2=631
        else:
            xmin1=82
            xmax1=313
            ymin=64
            ymax=298
            xmin2=360
            xmax2=594
        if(event.x>xmin1 and event.x<xmax1 and event.y>ymin and event.y<ymax):
            x=event.x-xmin1
            y=event.y-ymin
            i=int(math.floor(x*s.nz/(xmax1-xmin1+1)))
            j=int(math.floor(y*s.nz/(xmax1-xmin1+1)))
            if(i>=s.nz_min and i<=s.nz_max and j>=s.nz_min and j<=s.nz_max):
                i=i-s.nz_min
                j=j-s.nz_min
                if j!=i:
                    self.maskData[i,j]=not self.maskData[i,j]
                self.maskData[j,i]=not self.maskData[j,i]
                s.Rij_mask=np.copy(self.maskData)
                self.img=s.show_Rij(mask=True,normalization=True)
                s.Rij_mask=np.copy(self.original)
                self.drawCanvas()
        elif(event.x>xmin2 and event.x<xmax2 and event.y>ymin and event.y<ymax):
            x=event.x-xmin2
            y=event.y-ymin
            i=int(math.floor(x*s.nz/(xmax1-xmin1+1)))
            j=int(math.floor(y*s.nz/(xmax1-xmin1+1)))
            if(i>=s.nz_min and i<=s.nz_max and j>=s.nz_min and j<=s.nz_max):
                i=i-s.nz_min
                j=j-s.nz_min
                if j!=i:
                    self.maskData[i,j]=not self.maskData[i,j]
                self.maskData[j,i]=not self.maskData[j,i]
                s.Rij_mask=np.copy(self.maskData)
                self.img=s.show_Rij(mask=True,normalization=True)
                s.Rij_mask=np.copy(self.original)
                self.drawCanvas()
    def drawCanvas(self):
        canvas = FigureCanvasTkAgg(self.img,self.frame)
        canvas.show()
        canvas.get_tk_widget().grid(column=0,row=0)
        canvas.get_tk_widget().config(highlightthickness=0)
        canvas.get_tk_widget().bind("<Button-1>", self.OnMouseDown)

#Displays the plots 
def showView(root):
    currentx=0
    currenty=0
    global loaded
    if loaded:
        global canvas,frame
        figEh=plt.Figure()
        outerouterframe=tk.Frame()
        canvas=tk.Canvas(outerouterframe,width=int(10.23*figEh.get_dpi()),height=int(6.5*figEh.get_dpi()))
        ysb = tk.Scrollbar(outerouterframe,orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=ysb.set)
        ysb.pack(side="right",fill="y")
        frame=tk.Frame()
        canvas.pack(side="left", fill="both", expand=True)
        canvas.create_window(0,0,window=frame, anchor='nw')
        
        global averageImage
        figFrame=tk.Frame(frame)
        figFrame.grid(column=0,row=0,pady=5,padx=5)
        canvas1 = FigureCanvasTkAgg(averageImage, figFrame)
        tools=tkagg.NavigationToolbar2TkAgg(canvas1, figFrame)
        tools.update()
        canvas1.show()
        canvas1.get_tk_widget().config(highlightthickness=0)
        canvas1.get_tk_widget().pack()
        tools.pack()

        global fourierMask
        fourFrame=tk.Frame(frame)
        fourFrame.grid(column=1,row=0,pady=5,padx=5)
        canvas4=FigureCanvasTkAgg(fourierMask,fourFrame)
        tools2=tkagg.NavigationToolbar2TkAgg(canvas4,fourFrame)
        tools2.update()
        canvas4.show()
        canvas4.get_tk_widget().config(highlightthickness=0)
        canvas4.get_tk_widget().pack()
        tools2.pack()

        global rijMaskFalse
        canvas2 = FigureCanvasTkAgg(rijMaskFalse, master=frame)
        canvas2.get_tk_widget().config(highlightthickness=0)
        canvas2.show()
        canvas2.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
        canvas2.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)
        canvas2.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)
        
        global rijMaskTrue
        canvas3 = FigureCanvasTkAgg(rijMaskTrue, master=frame)
        canvas3.get_tk_widget().config(highlightthickness=0)
        canvas3.show()
        canvas3.get_tk_widget().grid(column=1, row=1,padx=5,pady=5)
        canvas3.get_tk_widget().bind("<ButtonPress-1>",rijTrueOnMouseDown)
        canvas3.get_tk_widget().bind("<ButtonRelease-1>",rijTrueOnMouseRelease)
        
        global sliceUn
        outerframe=tk.Frame(frame)
        imageFrame=tk.Frame(outerframe)
        backButton=tk.Button(outerframe, text="<", font=MEDIUMLARGE,command=lambda: backButtPress(imageFrame,s.imstack,"Stack Image "))
        forwardButton=tk.Button(outerframe,text=">",font=MEDIUMLARGE,command=lambda: forwardButtPress(imageFrame,s.imstack,"Stack Image "))
        backButton.pack(side="left",fill="y")
        forwardButton.pack(side="right",fill="y")
        imageFrame.pack(side="top",fill="both")
        outerframe.grid(row=2,column=0)
        forwardButtPress(imageFrame,s.imstack,"Stack Image ")
        
        global sliceReg
        outerframe2=tk.Frame(frame)
        imageFrame2=tk.Frame(outerframe2)
        backButton2=tk.Button(outerframe2, text="<", font=MEDIUMLARGE,command=lambda: backButtPress(imageFrame2,s.stack_registered,"Registered Image "))
        forwardButton2=tk.Button(outerframe2,text=">",font=MEDIUMLARGE,command=lambda: forwardButtPress(imageFrame2,s.stack_registered,"Registered Image "))
        backButton2.pack(side="left",fill="y")
        forwardButton2.pack(side="right",fill="y")
        imageFrame2.pack(side="top",fill="both")
        outerframe2.grid(row=2,column=1)
        forwardButtPress(imageFrame2,s.stack_registered,"Registered Image ")
        
        root.update()
        outerouterframe.grid(column=0,row=0)
        canvas.config(scrollregion=canvas.bbox("all"))

#For moving forward in the slice viewer        
def forwardButtPress(frame,stack,txt):
    if txt=="Stack Image ":
        global sliceUn
        if sliceUn==s.nz-1:
            sliceUn=0
        else: 
            sliceUn=sliceUn+1
        slice=sliceUn
    else:
        global sliceReg
        if sliceReg==s.nz-1:
            sliceReg=0
        else:
            sliceReg=sliceReg+1
        slice=sliceReg
    fig,ax=plt.subplots(figsize=(4,4),dpi=100)
    ax.matshow(stack[:,:,slice],cmap='gray')
    ax.axis('off')
    ax.set_title(txt+str(slice))
    can=FigureCanvasTkAgg(fig,frame)
    can.get_tk_widget().config(highlightthickness=0)
    can.get_tk_widget().grid(column=0,row=0)
    return slice
    
#For moving backward in the slice viewer
def backButtPress(frame,stack,txt):
    if txt=="Stack Image ":
        global sliceUn
        if sliceUn==0:
            sliceUn=s.nz-1
        else: 
            sliceUn=sliceUn-1
        slice=sliceUn
    else:
        global sliceReg
        if sliceReg==0:
            sliceReg=s.nz-1
        else:
            sliceReg=sliceReg-1
        slice=sliceReg
    fig,ax=plt.subplots(figsize=(4,4),dpi=100)
    ax.matshow(stack[:,:,slice],cmap='gray')
    ax.axis('off')
    ax.set_title(txt+str(slice))
    can=FigureCanvasTkAgg(fig,frame)
    can.get_tk_widget().config(highlightthickness=0)
    can.get_tk_widget().grid(column=0,row=0)
    return slice
    
#Mask=False Rij Mouse Press
def rijFalseOnMouseDown(event):
    global frame
    canvas2 = FigureCanvasTkAgg(s.show_Rij(mask=True,normalization=False),master=frame)
    canvas2.get_tk_widget().config(highlightthickness=0)
    canvas2.show()
    canvas2.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
    canvas2.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)
    canvas2.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)

#Mask=False Rij Mouse Release
def rijFalseOnMouseRelease(event):
    global rijMaskFalse,frame
    canvas2 = FigureCanvasTkAgg(rijMaskFalse,master=frame)
    canvas2.get_tk_widget().config(highlightthickness=0)
    canvas2.show()
    canvas2.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
    canvas2.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)
    canvas2.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)

#Mask=True Rij Mouse Press
def rijTrueOnMouseDown(event):
    global frame
    canvas3 = FigureCanvasTkAgg(s.show_Rij(mask=True,normalization=False), master=frame)
    canvas3.get_tk_widget().config(highlightthickness=0)
    canvas3.show()
    canvas3.get_tk_widget().grid(column=1, row=1,padx=5,pady=5)
    canvas3.get_tk_widget().bind("<ButtonRelease-1>",rijTrueOnMouseRelease)
    canvas3.get_tk_widget().bind("<ButtonPress-1>",rijTrueOnMouseDown)

#Mask=True Rij Mouse Release
def rijTrueOnMouseRelease(event):
    global frame
    canvas3 = FigureCanvasTkAgg(rijMaskTrue, master=frame)
    canvas3.get_tk_widget().config(highlightthickness=0)
    canvas3.show()
    canvas3.get_tk_widget().grid(column=1, row=1,padx=5,pady=5)
    canvas3.get_tk_widget().bind("<ButtonRelease-1>",rijTrueOnMouseRelease)
    canvas3.get_tk_widget().bind("<ButtonPress-1>",rijTrueOnMouseDown)
    
#Recalculates data without redoing entire thing if parameter is changed. Saves time. Good thing
def recalc():
    s.get_averaged_image(crop=False)
    global rijMaskFalse
    rijMaskFalse = s.show_Rij(mask=False,normalization=False)
    global rijMaskTrue
    rijMaskTrue = s.show_Rij(mask=True,normalization=True)
    global averageImage,fourierMask
    fourierMask=s.show_Fourier_mask(image_index=0)
    averageImage = s.show()
    showView(root)



#Runs s.get_outliers with correct arguments        
def outliers():
    global outlierMethod
    if outlierMethod=="NN":
        global shiftMax
        s.get_outliers("NN",shiftMax)
    elif outlierMethod=="PF":
        global polycoeffmax, zscore
        s.get_outliers("PF",polycoeffmax,zscore)
    
#Runs when you select your data set
def loadData(fileType):
    t0=time()
    global s
    if fileType=='tif':
        filepath=os.path.dirname(askopenfilename())
        files=[]
        if(fileType=='tif'):
            for filename in listdir(filepath):
                files.append(filepath+"/"+filename)
        im=np.array(Image.open(files[0]))
        fov=im.shape[0]
        stack=np.empty((fov,fov,len(files)))
        for i in range(len(files)):
            im=np.array(Image.open(files[i]))
            im=im/float(2**16)
            stack[:,:,i]=im	
    elif fileType=='ser':
        stack=serReader.serReader(askopenfilename())
    else:
        stack=np.rollaxis(imread(askopenfilename()),0,3)
        stack=stack[:,:,:]/float(2**16)
        
    # Instantiate imstack object
    s =stackregister.imstack(stack[:,:,:])	 
    global correlationType,findMaxima,loaded,fouriern,fourierMaskType,canvas
    #canvas=None
    s.getFFTs()
    s.makeFourierMask(mask=fourierMaskType,n=fouriern)
    s.findImageShifts(correlationType=correlationType,findMaxima=findMaxima,verbose=False)
    t=time()-t0
    s.set_nz(0,s.nz)
    loaded=True
    outliers()
    recalc()
    
    print("Completed calculations in {} minutes {} seconds".format(int(t/60),t%60))
   
        
#Saves the shift matrixes to a text file
def saveTxt():
    fileName=tkFileDialog.asksaveasfilename(defaultextension='.txt')
    with file(fileName, 'w') as outfile:
        outfile.write('X shifts\n')
        np.savetxt(outfile, s.X_ij, fmt='%-7.2f')
        outfile.write('\nY Shifts\n')
        np.savetxt(outfile, s.Y_ij, fmt='%-7.2f')
    global loaded
    if loaded == False:
        tkMessageBox.showwarning("Load Data","Please load images in to analze first")

#This redoes the necessary parts if you change the outlier method    
def outRedo(bool,coef,z,shift):
    print bool
    global polycoeffmax,zscore,shiftMax,outlierMethod
    polycoeffmax = int(coef)
    zscore = float(z)
    shiftMax = float(shift)
    if not bool:
        outlierMethod = "NN"
    else:
        outlierMethod="PF"
    outliers()
    global outpop
    outpop.destroy()
    recalc()
    
#This redoes the necessary parts if you change the nz range
def nzRedo(minText,maxText,listEntry):
    s.set_nz(int(minText),int(maxText))
    s.set_bad_images(map(int,listEntry.split(",")))
    global nzpop
    nzpop.destroy()
    outliers()
    recalc()

#Updates fourier mask and reruns all necessary calculations after changing it    
def fourierRedo(n,type):
    global fouriern, fourierpop,fourierMaskType
    fouriern=int(n)
    fourierMaskType=type
    s.makeFourierMask(mask=type,n=fouriern)
    s.findImageShifts(correlationType=correlationType,findMaxima=findMaxima,verbose=False)
    outliers()
    recalc()
 
#Updates the method used to find shifts between images and reruns all necessary calculations 
def shiftUpdate(type,maxima,numPeak,sigma,windowRad,numiter,minWinFrac):
    global shiftpop, correlationType, findMaxima, gaussiannumpeaks,sigmaguess,windowradius,minwindowfrac
    correlationType=type
    findMaxima=maxima
    gaussiannumpeaks=numPeak
    sigmaguess=sigma
    windowradius=windowRad
    numiter=numiter
    minwindowfrac=minWinFrac
    s.setGaussianFitParams(num_peaks=numPeak,sigma_guess=sigma,window_radius=windowRad)
    s.setCoMParams(num_iter=numiter,min_window_frac=minWinFrac)
    s.findImageShifts(correlationType=type,findMaxima=maxima,verbose=False)
    shiftpop.destroy()
    recalc()
 
#This makes the radio buttons in the change outlier method work
def radioUpdate(bool,outCoeEntry,outZEntry,outShiftEntry):
    if bool:
        outCoeEntry.configure(state=tk.NORMAL)
        outZEntry.configure(state=tk.NORMAL)
        outShiftEntry.configure(state=tk.DISABLED)
    else:
        outCoeEntry.configure(state=tk.DISABLED)
        outZEntry.configure(state=tk.DISABLED)
        outShiftEntry.configure(state=tk.NORMAL)

#This makes the radio buttons in the Change Shift Method work        
def otherRadioUpdate(selection,numPeakEntry,sigmaEntry,windowRadEntry,numiterEntry,minWinFracEntry):
    numPeakEntry.configure(state=tk.DISABLED)
    sigmaEntry.configure(state=tk.DISABLED)
    windowRadEntry.configure(state=tk.DISABLED)
    numiterEntry.configure(state=tk.DISABLED)
    minWinFracEntry.configure(state=tk.DISABLED)
    if selection=="gf":
        numPeakEntry.configure(state=tk.NORMAL)
        sigmaEntry.configure(state=tk.NORMAL)
        windowRadEntry.configure(state=tk.NORMAL)
    elif selection=="com":
        numiterEntry.configure(state=tk.NORMAL)
        minWinFracEntry.configure(state=tk.NORMAL)
               
#This creates the popup to change the settings when you change outlier method
def outlierPopup():
    global loaded
    if loaded:
        global outpop
        outpop = tk.Tk()
        outpop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
        outframe = ttk.Frame(outpop, padding="3 3 12 12")
        outframe.grid(column=0,row=0)
        outframe.columnconfigure(0,weight=1)
        outframe.rowconfigure(0,weight=1)
        outpop.wm_title("Edit Outlier Criteria")
        
        v=tk.BooleanVar(outpop)
        global outlierMethod
        if outlierMethod=="PF":
            v.set(True)
        else:
            v.set(False)
        
        fitframe = tk.LabelFrame(outframe,text='Polynomial Fit Settings',relief='sunken')
        aroundframe = tk.LabelFrame(outframe,text='Nearest Points Fit Settings',relief='sunken')
        fitframe.grid(column=0,row=1,columnspan=2)
        aroundframe.grid(column=0,row=2,columnspan=2)
        
        global polycoeffmax,zscore
        outCoeLabel=tk.Label(fitframe, text="Coefficent of Polynomial fit:",font=LARGE_FONT)
        outCoeLabel.grid(column=0,row=0)
        outCoeEntry = tk.Entry(fitframe, width=7, font=LARGE_FONT)
        outCoeEntry.grid(column=1, row=0)
        outCoeEntry.insert(0,polycoeffmax)
        outZLabel = tk.Label(fitframe, text="Z Value:",font=LARGE_FONT)
        outZLabel.grid(column=0, row=1)
        outZEntry = tk.Entry(fitframe, width=7, font=LARGE_FONT)
        outZEntry.grid(column=1, row=1)
        outZEntry.insert(0,zscore)
        
        global shiftMax
        outShiftLabel=tk.Label(aroundframe,text="Max Shift:",font=LARGE_FONT)
        outShiftLabel.grid(column=0,row=0)
        outShiftEntry = tk.Entry(aroundframe,width=7,font=LARGE_FONT)
        outShiftEntry.grid(column=1,row=0)
        outShiftEntry.insert(0,shiftMax)
        
        outFitRadio=tk.Radiobutton(outframe, text="Use Polynomial Fit", variable=v, value=True,command=lambda: radioUpdate(True,outCoeEntry,outZEntry,outShiftEntry))
        outAroundRadio=tk.Radiobutton(outframe, text="Use nearest points", variable=v, value=False,command=lambda: radioUpdate(False,outCoeEntry,outZEntry,outShiftEntry))
        radioUpdate(v.get(),outCoeEntry,outZEntry,outShiftEntry)
        outFitRadio.grid(column=0, row=0)
        outAroundRadio.grid(column=1,row=0)
        if v.get():
            outFitRadio.select()
        else:
            outAroundRadio.select()
        outSaveButton = tk.Button(outframe, text="Save", font=MEDIUMLARGE,command=lambda: outRedo(v.get(),outCoeEntry.get(),outZEntry.get(),outShiftEntry.get()))
        outSaveButton.grid(column=0,row=3)
        outCancelButton = tk.Button(outframe, text="Cancel", font=MEDIUMLARGE, command=outpop.destroy)
        outCancelButton.grid(column=1,row=3)
        
        outpop.mainloop()
    else:
       tkMessageBox.showwarning("Load Data","Please load images in to analyze first.") 
       
#This makes the popup to change the settings for the nz range
def nzPopup():
    global loaded
    if loaded:
        global nzpop
        nzpop = tk.Tk()
        nzpop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
        nzframe = ttk.Frame(nzpop, padding="3 3 12 12")
        nzframe.grid(column=0,row=0)
        nzframe.columnconfigure(0,weight=1)
        nzframe.rowconfigure(0,weight=1)
        nzpop.wm_title("Edit nz range")
        nzminLabel=tk.Label(nzframe, text="nz min:",font=LARGE_FONT)
        nzminLabel.grid(column=0,row=0)
        nzminEntry = tk.Entry(nzframe, width=7, font=LARGE_FONT)
        nzminEntry.grid(column=1, row=0)
        nzminEntry.insert(0,str(s.nz_min))
        nzmaxLabel = tk.Label(nzframe, text="nz max:",font=LARGE_FONT)
        nzmaxLabel.grid(column=0, row=1)
        nzmaxEntry = tk.Entry(nzframe, width=7, font=LARGE_FONT)
        nzmaxEntry.grid(column=1, row=1)
        nzmaxEntry.insert(0,str(s.nz_max))
        nzListLabel = tk.Label(nzframe,text="Bad images:",font=LARGE_FONT)
        nzListLabel.grid(column=0,row=2)
        nzListEntry = tk.Entry(nzframe,width=7,font=LARGE_FONT)
        nzListEntry.grid(column=1,row=2)
        if hasattr(s, 'bad_images'):
            nzListEntry.insert(0,','.join(map(str, s.bad_images)))
        nzSaveButton = tk.Button(nzframe, text="Save", font=MEDIUMLARGE,command=lambda: nzRedo(nzminEntry.get(),nzmaxEntry.get(),nzListEntry.get()))
        nzSaveButton.grid(column=0,row=3)
        nzCancelButton = tk.Button(nzframe, text="Cancel", font=MEDIUMLARGE, command=nzpop.destroy)
        nzCancelButton.grid(column=1,row=3)
    else:
        tkMessageBox.showwarning("Load Data","Please load images in to analze first")

#This makes the popup to change if there is a fourier transform mask or not
def fourierPopup():
    global loaded
    if loaded:
        global fourierpop, fouriern,fourierMaskType
        fourierpop= tk.Tk()
        fourierpop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
        fourierpop.title("Edit Fourier Transform Mark")
        fourierNLabel = tk.Label(fourierpop,text="N Value:",font=LARGE_FONT)
        fourierNLabel.grid(column=0,row=0)
        fourierNEntry=tk.Entry(fourierpop,width=2,font=LARGE_FONT)
        fourierNEntry.grid(column=1,row=0)
        fourierNEntry.insert(0,fouriern)
        
        radioText = tk.StringVar(fourierpop)
        radioText.set(fourierMaskType)
        fourierRadioFrame=tk.LabelFrame(fourierpop,text="Select Mask Type", relief="sunken")
        fourierNoneRadio=tk.Radiobutton(fourierRadioFrame,text="None",variable=radioText,value="none")
        fourierLowRadio=tk.Radiobutton(fourierRadioFrame,text="Lowpass",variable=radioText,value="lowpass")
        fourierBandRadio=tk.Radiobutton(fourierRadioFrame,text="Bandpass",variable=radioText,value="bandpass")
        fourierRadioFrame.grid(column=0,columnspan=2,row=1)
        fourierNoneRadio.grid(column=0,row=0)
        fourierLowRadio.grid(column=1,row=0)
        fourierBandRadio.grid(column=2,row=0)
        if radioText.get()=="none":
            fourierNoneRadio.select()
        elif radioText.get()=="lowpass":
            fourierLowRadio.select()
        else:
            fourierBandRadio.select()
        fourierSaveButton=tk.Button(fourierpop,text="Save",font=MEDIUMLARGE,command=lambda: fourierRedo(fourierNEntry.get(),radioText.get()))
        fourierSaveButton.grid(column=0,row=2)
        fourierCancelButton=tk.Button(fourierpop,text="Cancel",font=MEDIUMLARGE,command=fourierpop.destroy)
        fourierCancelButton.grid(column=1,row=2)
        fourierpop.mainloop()
    else:
        tkMessageBox.showwarning("Load Data","Please load images in to analyze first.")

#This one makes a popup to change the way it calculates the shifts between the images (yay!)
def shiftPopup():
    global shiftpop, correlationType,findMaxima,gaussiannumpeaks,sigmaguess,windowradius,numiter,minwindowfrac
    shiftpop=tk.Tk()
    shiftpop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
    shiftpop.title("Edit Image Shift Calculation")    
    typeFrame=tk.LabelFrame(shiftpop,text="Select Correleation Type",relief="ridge")
    typeRadio=tk.StringVar(shiftpop)
    typeRadio.set(correlationType)
    typeFrame.grid(row=0,column=0,columnspan=2)
    shiftCrossRadio=tk.Radiobutton(typeFrame,text="Cross Correleation",variable=typeRadio,value="cc")
    shiftMutualRadio=tk.Radiobutton(typeFrame,text="Mutual Correlation",variable=typeRadio,value="mc")
    shiftPhaseRadio=tk.Radiobutton(typeFrame,text="Phase Correleation",variable=typeRadio,value="pc")
    shiftCrossRadio.grid(column=0,row=0)
    shiftMutualRadio.grid(column=1,row=0)
    shiftPhaseRadio.grid(column=2,row=0)
    if correlationType=="cc":
        shiftCrossRadio.select()
    elif correlationType=="mc":
        shiftMutualRadio.select()
    else:
        shiftPhaseRadio.select()
    
    maximaFrame=tk.LabelFrame(shiftpop,text="Select Method for Finding Maxima",relief="ridge")
    maximaFrame.grid(row=1,column=0,columnspan=2)
    maximaRadio=tk.StringVar(shiftpop)
    maximaRadio.set(findMaxima)
    maximaPixelRadio=tk.Radiobutton(maximaFrame,text="Pixel",variable=maximaRadio,value="pixel",command=lambda: otherRadioUpdate("pixel",numPeakEntry,sigmaEntry,windowRadEntry,numiterEntry,minWinFracEntry))
    maximaGfRadio=tk.Radiobutton(maximaFrame,text="Subpixel by Gaussian Fit",variable=maximaRadio,value="gf",command=lambda: otherRadioUpdate("gf",numPeakEntry,sigmaEntry,windowRadEntry,numiterEntry,minWinFracEntry))
    #maximaComRadio=tk.Radiobutton(maximaFrame,text="Center of Mass",variable=maximaRadio,value="com",command=lambda: otherRadioUpdate("com",numPeakEntry,sigmaEntry,windowRadEntry,numiterEntry,minWinFracEntry))
    maximaPixelRadio.grid(column=0,row=0)
    maximaGfRadio.grid(column=1,row=0)
    #maximaComRadio.grid(column=2,row=0)
    if findMaxima=="pixel":
        maximaPixelRadio.select()
    elif findMaxima=="gf":
        maximaGfRadio.select()
    else:
        maximaComRadio.select()    
    
    gfSettingFrame=tk.LabelFrame(shiftpop,text="Gaussian Fit Settings",relief="sunken")
    gfSettingFrame.grid(column=0,row=2)
    numPeakLabel=tk.Label(gfSettingFrame,text="Number of Peaks:",font=MEDIUMLARGE)
    numPeakEntry=tk.Entry(gfSettingFrame,width=2,font=LARGE_FONT)
    numPeakEntry.insert(0,str(gaussiannumpeaks))
    sigmaLabel=tk.Label(gfSettingFrame,text="Sigma Guess:",font=MEDIUMLARGE)
    sigmaEntry=tk.Entry(gfSettingFrame,width=2,font=LARGE_FONT)
    sigmaEntry.insert(0,str(sigmaguess))
    windowRadLabel=tk.Label(gfSettingFrame,text="Window Radius:",font=MEDIUMLARGE)
    windowRadEntry=tk.Entry(gfSettingFrame,width=2,font=LARGE_FONT)
    windowRadEntry.insert(0,str(windowradius))
    numPeakLabel.grid(column=0,row=0)
    numPeakEntry.grid(column=1,row=0)
    sigmaLabel.grid(column=0,row=1)
    sigmaEntry.grid(column=1,row=1)
    windowRadLabel.grid(column=0,row=2)
    windowRadEntry.grid(column=1,row=2)
    
    comSettingFrame=tk.LabelFrame(shiftpop,text="Center of Mass Settings",relief="sunken")
    comSettingFrame.grid(column=1,row=2)
    numiterLabel=tk.Label(comSettingFrame,text="Number of Iterations:",font=MEDIUMLARGE)
    numiterEntry=tk.Entry(comSettingFrame,width=2,font=LARGE_FONT)
    numiterEntry.insert(0,str(numiter))
    minWinFracLabel=tk.Label(comSettingFrame,text="Minimum Window Fraction:",font=MEDIUMLARGE)
    minWinFracEntry=tk.Entry(comSettingFrame,width=2,font=LARGE_FONT)
    minWinFracEntry.insert(0,str(minwindowfrac))
    numiterLabel.grid(column=0,row=0)
    numiterEntry.grid(column=1,row=0)
    minWinFracLabel.grid(column=0,row=1)
    minWinFracEntry.grid(column=1,row=1)
    
    shiftSaveButton=tk.Button(shiftpop,text="Save",font=MEDIUMLARGE,command=lambda: shiftUpdate(typeRadio.get(),maximaRadio.get(),numPeakEntry.get(),sigmaEntry.get(),windowRadEntry.get(),numiterEntry.get(),minWinFracEntry.get()))
    shiftSaveButton.grid(column=0,row=3)
    shiftCancelButton=tk.Button(shiftpop,text="Cancel",font=MEDIUMLARGE,command=shiftpop.destroy)
    shiftCancelButton.grid(column=1,row=3)
    otherRadioUpdate(maximaRadio.get(),numPeakEntry,sigmaEntry,windowRadEntry,numiterEntry,minWinFracEntry)
    shiftpop.mainloop()
    
#This makes the great about popup
def aboutPopup():
    aboutPop=tk.Tk()
    aboutPop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
    aboutPop.title("About")
    aboutPop.resizable(width=False, height=False)
    t1=tk.Label(aboutPop,text="Grr (Great Rigid Registration)")
    t1.grid(column=0,row=0)
    t2=tk.Label(aboutPop,text="v1.0   Released # Month 2017")
    t2.grid(column=0,row=1)
    t3=tk.Label(aboutPop,text="")
    t3.grid(column=0,row=2)
    t4=tk.Label(aboutPop,text="Emily Waite, Ismail El Baggari, Benjamin Savitzky, Lena F. Kourkoutis")
    t4.grid(column=0,row=3)
    t5=tk.Label(aboutPop,text="Department of Applied and Engineering Physics, Cornell University")
    t5.grid(column=0,row=4)

#Opens user Guide pdf 
def openHelp():
    filename=os.path.join(os.path.dirname(os.path.abspath(__file__)),"User Guide.pdf")
    if sys.platform == "win32":
        os.startfile(filename)
    else:
        opener ="open" if sys.platform == "darwin" else "xdg-open"
        subprocess.call([opener, filename]) 

#Makes the basic window and the menu bars you need
def rootSetup(root):
    root.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
    root.title('Grr')
    fig=plt.Figure()
    root.maxsize(width=int(11.5*fig.get_dpi()),height=int(7.5*fig.get_dpi()))
    root.resizable(width=False, height=False)
    menubar = tk.Menu(root)
    filemenu=tk.Menu(menubar, tearoff=0)
    filemenu.add_command(label="Import Extracted .tif files", command=lambda: loadData(fileType='tif'))
    filemenu.add_command(label="Import .tif stack",command=lambda:loadData(fileType='tifstack'))
    filemenu.add_command(label="Import .ser File",command=lambda: loadData(fileType='ser')) 
    filemenu.add_command(label="Save Image")
    filemenu.add_command(label="Generate Report")
    filemenu.add_command(label="Export Shift Data")
    filemenu.entryconfig(3,command=lambda: s.save(tkFileDialog.asksaveasfilename(defaultextension='.tif')))
    filemenu.entryconfig(4, command=lambda: s.save_report(tkFileDialog.asksaveasfilename(defaultextension='.pdf')))
    filemenu.entryconfig(5,command=saveTxt)
    editmenu=tk.Menu(menubar, tearoff=0)
    editmenu.add_command(label="Change nz range",command=nzPopup)
    editmenu.add_command(label="Change Outlier Method", command=outlierPopup)
    editmenu.add_command(label="Custom mask",command=lambda: CustomMask(master=root))
    editmenu.add_command(label="Fourier Transformation Mask",command=fourierPopup)
    editmenu.add_command(label="Change Image Shift Method",command=shiftPopup) 
    helpmenu=tk.Menu(menubar, tearoff=0)
    helpmenu.add_command(label="About",command=aboutPopup)
    helpmenu.add_command(label="Getting Started",command=openHelp) 
    menubar.add_cascade(label="File", menu=filemenu)
    menubar.add_cascade(label="Edit", menu=editmenu)
    menubar.add_cascade(label="Help", menu=helpmenu)
    root.config(menu=menubar)        
    
root = tk.Tk()
rootSetup(root)
root.mainloop()  