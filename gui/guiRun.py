#!/usr/bin/env python

"""
guiRun.py
A Graphical User Interphase to be used with rigidRegistration library

"""

#Import global libraries
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.backends.backend_tkagg as tkagg
from time import time
from PIL import Image
from os import listdir
import subprocess
from tifffile import imread
import math as math
import sys
import os
if sys.version_info[0] < 3:
    import Tkinter as tk
    import Tkinter.ttk as ttk
    from Tkinter.tkFileDialog import askopenfilename #ttk overwrites this one 
    import Tkinter.tkFileDialog
    import Tkinter.Tkconstants
    import Tkinter.tkMessageBox

else:
    import tkinter as tk
    import tkinter.ttk as ttk
    from tkinter.filedialog import askopenfilename #ttk overwrites this one 
    import tkinter.filedialog as tkFileDialog
    import tkinter.constants as Tkconstats
    import tkinter.messagebox as tkMessageBox

#Finds and imports local libraries    
sys.path.append('../')
from rigidregistration import stackregistration 
import serReader

#Initialize global variables (will be given values later)
s=None
rijMaskFalse=None
rijMaskTrue=None
averageImage=None
fourierMask=None
mainCanvas=None
frame=None
outpop = None
nzpop = None
fourierpop = None
outCoeEntry = None
outZEntry = None
outShiftEntry = None
imageFrameReg=None

loaded=False
sliceUn=-1
sliceReg=-1
stackUsed=None
MEDIUMLARGE = ("Verdana",20)
MEDIUMFONT= ("Verdana",14)
SMALLFONT=("Verdana",12)

########## Methods for building actual Interface ##########

class CustomToolbar(tkagg.NavigationToolbar2TkAgg):
    """
    Edits the toolbar library that is used on some of the figures to fit our purpose
    """
    def __init__(self,canvas_,parent_):
        self.toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Back', 'Back to  previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            (None, None, None, None),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
            (None, None, None, None),
            (None, None, None, None),
            ('Save', 'Save the figure', 'filesave', 'save_figure'),
            )
        tkagg.NavigationToolbar2TkAgg.__init__(self,canvas_,parent_)

def loadParams():
    """
    Loads default values for the calculations from defaultParam.txt
    """
    filename=os.path.join(os.path.dirname(os.path.abspath(__file__)),"defaultParam.txt")
    global threshold, maxpaths, findMaxima,fouriern,fourierMaskType,gaussiannumpeaks,sigmaguess,windowradius
    if os.path.isfile(filename):
        import imp
        f = open(filename)
        data = imp.load_source('data', '', f)
        f.close()
        threshold=data.threshold
        maxpaths=data.maxpaths
        findMaxima=data.findMaxima
        fouriern=data.fouriern
        fourierMaskType=data.fourierMaskType
        gaussiannumpeaks=data.gaussiannumpeaks
        sigmaguess=data.sigmaguess
        windowradius=data.windowradius
        
    else:
        threshold=10
        maxpaths=5
        findMaxima="pixel"
        fouriern = 4
        fourierMaskType="bandpass"
        gaussiannumpeaks=3
        sigmaguess=2
        windowradius=6      
         
def showView(root):
    """
    Puts all the plots correctly into the window when all results are ready
    """
    currentx=0
    currenty=0
    global loaded
    if loaded:
        global mainCanvas,frame
        figOut=plt.Figure()
        outerouterframe=tk.Frame()
        mainCanvas=tk.Canvas(outerouterframe,width=int(10.23*figOut.get_dpi()),height=int(6.5*figOut.get_dpi()))
        ysb = tk.Scrollbar(outerouterframe,orient="vertical", command=mainCanvas.yview)
        mainCanvas.configure(yscrollcommand=ysb.set)
        ysb.pack(side="right",fill="y")
        frame=tk.Frame()
        mainCanvas.pack(side="left", fill="both", expand=True)
        mainCanvas.create_window(0,0,window=frame, anchor='nw')
        
        global averageImage
        figFrame=tk.Frame(frame)
        figFrame.grid(column=0,row=0,pady=5,padx=5)
        avgImgCanvas = FigureCanvasTkAgg(averageImage, figFrame)
        imgTools=CustomToolbar(avgImgCanvas, figFrame)
        imgTools.update()
        avgImgCanvas.show()
        avgImgCanvas.get_tk_widget().config(highlightthickness=0)
        avgImgCanvas.get_tk_widget().pack()
        imgTools.pack()

        global fourierMask
        fourFrame=tk.Frame(frame)
        fourFrame.grid(column=1,row=0,pady=5,padx=5)
        fourierCanvas=FigureCanvasTkAgg(fourierMask,fourFrame)
        fourierTools=CustomToolbar(fourierCanvas,fourFrame)
        fourierTools.update()
        fourierCanvas.show()
        fourierCanvas.get_tk_widget().config(highlightthickness=0)
        fourierCanvas.get_tk_widget().pack()
        fourierTools.pack()

        global rijMaskFalse
        rijFalseCanvas = FigureCanvasTkAgg(rijMaskFalse, master=frame)
        rijFalseCanvas.get_tk_widget().config(highlightthickness=0)
        rijFalseCanvas.show()
        rijFalseCanvas.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
        rijFalseCanvas.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)
        rijFalseCanvas.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)
        
        global rijMaskTrue
        rijTrueCanvas = FigureCanvasTkAgg(rijMaskTrue, master=frame)
        rijTrueCanvas.get_tk_widget().config(highlightthickness=0)
        rijTrueCanvas.show()
        rijTrueCanvas.get_tk_widget().grid(column=1, row=1,padx=5,pady=5)
        rijTrueCanvas.mpl_connect('button_press_event',rijTrueOnMouseDown)
        
        global sliceUn
        outerframe=tk.Frame(frame)
        imageFrame=tk.Frame(outerframe)
        backButton=tk.Button(outerframe, text="<", font=MEDIUMLARGE,command=lambda: backButtPress(imageFrame,s.imstack,"Stack Image ",lng=s.nz))
        forwardButton=tk.Button(outerframe,text=">",font=MEDIUMLARGE,command=lambda: forwardButtPress(imageFrame,s.imstack,"Stack Image ",lng=s.nz))
        backButton.pack(side="left",fill="y")
        forwardButton.pack(side="right",fill="y")
        imageFrame.pack(side="top",fill="both")
        outerframe.grid(row=2,column=0)
        forwardButtPress(imageFrame,s.imstack,"Stack Image ")
        
        global sliceReg,stackUsed,imageFrameReg,forwardButton2,backButton2
        outerframe2=tk.Frame(frame)
        imageFrameReg=tk.Frame(outerframe2)
        stackUsed=s.stack_registered
        backButton2=tk.Button(outerframe2, text="<", font=MEDIUMLARGE,command=lambda: backButtPress(imageFrameReg,stackUsed,"Registered Image ",lng=s.nz))
        forwardButton2=tk.Button(outerframe2,text=">",font=MEDIUMLARGE,command=lambda: forwardButtPress(imageFrameReg,stackUsed,"Registered Image ",lng=s.nz))
        backButton2.pack(side="left",fill="y")
        forwardButton2.pack(side="right",fill="y")
        imageFrameReg.pack(side="top",fill="both")
        outerframe2.grid(row=2,column=1)
        forwardButtPress(imageFrameReg,stackUsed,"Registered Image ")
        
        root.update()
        outerouterframe.grid(column=0,row=0)
        mainCanvas.config(scrollregion=mainCanvas.bbox("all"))
        
        root.bind("<Button-4>", mouse_wheel)
        root.bind("<Button-5>", mouse_wheel)
        root.bind("<MouseWheel>", mouse_wheel)

def mouse_wheel(event):
    """
    Binds the mouse scrollwheel to the vertical scrollbar
    """
    global mainCanvas
    if sys.platform == "darwin":
        mainCanvas.yview_scroll(int(-1*(event.delta)), "units")
    else:
        mainCanvas.yview_scroll(int(-1*(event.delta/120)), "units")
 
######## Methods for interacting with buttons in Interface ##########
 
def forwardButtPress(frame,stack,txt,lng=0):
    """
    This moves forward to the next slice in the slice viewer
    when the forward button is pressed
    """
    if txt=="Stack Image ":
        global sliceUn
        if sliceUn==lng-1:
            sliceUn=0
        else: 
            sliceUn=sliceUn+1
        slice=sliceUn
    else:
        global sliceReg
        if sliceReg==lng-1:
            sliceReg=0
        else:
            sliceReg=sliceReg+1
        slice=sliceReg
    fig,ax=plt.subplots(figsize=(4,4),dpi=100)
    ax.matshow(stack[:,:,int(slice)],cmap='gray')
    ax.axis('off')
    ax.set_title(txt+str(slice),y=1)
    if lng==2:
        ax.set_title(txt,y=1)
    plt.subplots_adjust(left=0,right=1,bottom=0.01,top=0.93,wspace=0.03,hspace=0.01)
    can=FigureCanvasTkAgg(fig,frame)
    can.get_tk_widget().config(highlightthickness=0)
    can.get_tk_widget().grid(column=0,row=0)
    return slice
    
def backButtPress(frame,stack,txt,lng=0):
    """
    This moves backwards to the previous slice in the slice viewer
    when the back button is pressed
    """
    if txt=="Stack Image ":
        global sliceUn
        if sliceUn==0:
            sliceUn=lng-1
        else: 
            sliceUn=sliceUn-1
        slice=sliceUn
    else:
        global sliceReg
        if sliceReg==0:
            sliceReg=lng-1
        else:
            sliceReg=sliceReg-1
        slice=sliceReg
    fig,ax=plt.subplots(figsize=(4,4),dpi=100)
    ax.matshow(stack[:,:,slice],cmap='gray')
    ax.axis('off')
    ax.set_title(txt+str(slice),y=1)
    if lng==2:
        ax.set_title(txt,y=1)
    plt.subplots_adjust(left=0,right=1,bottom=0.01,top=0.93,wspace=0.03,hspace=0.01)
    can=FigureCanvasTkAgg(fig,frame)
    can.get_tk_widget().config(highlightthickness=0)
    can.get_tk_widget().grid(column=0,row=0)
    return slice
    
def rijFalseOnMouseDown(event):
    """
    Defines the event when the mouse is clicked on the unmasked rij matrixes
    This will replace the unmasked rij with a normalized masked rij
    """
    global frame
    rijFalseCanvas = FigureCanvasTkAgg(s.show_Rij(mask=True,normalization=False,returnfig=True),master=frame)
    rijFalseCanvas.get_tk_widget().config(highlightthickness=0)
    rijFalseCanvas.show()
    rijFalseCanvas.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
    rijFalseCanvas.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)
    rijFalseCanvas.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)

def rijFalseOnMouseRelease(event):
    """
    Defines event when mouse button is released on originally unmasked rij
    Converts the normalized masked rij back to original unmasked rij
    """
    global rijMaskFalse,frame
    rijFalseCanvas = FigureCanvasTkAgg(rijMaskFalse,master=frame)
    rijFalseCanvas.get_tk_widget().config(highlightthickness=0)
    rijFalseCanvas.show()
    rijFalseCanvas.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
    rijFalseCanvas.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)
    rijFalseCanvas.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)

def rijTrueOnMouseDown(event):
    """
    Defines event for when mouse is presed down on the masked rij
    This will select whatever 2 slices are at the point where the user selected
    and display just these 2 slices in the slice viewer underneath
    """
    global stackUsed,sliceReg,imageFrameReg,backButton2,forwardButton2
    if event.inaxes is not None:
        i=int(round(event.xdata))
        j=int(round(event.ydata))
        sliceReg=0
        stackUsed=np.dstack((s.stack_registered[:,:,i],s.stack_registered[:,:,j]))
        strng="Comparing "+str(i)+" and "+str(j)
        backButton2.configure(command=lambda:backButtPress(imageFrameReg,stackUsed,strng,lng=2))
        forwardButton2.configure(command=lambda:forwardButtPress(imageFrameReg,stackUsed,strng,lng=2))
        forwardButtPress(imageFrameReg,stackUsed,strng,lng=2)
    else:
        forwardButton2.configure(command=lambda:forwardButtPress(imageFrameReg,stackUsed,"Registered Image ",lng=s.nz))
        backButton2.configure(command=lambda:backButtPress(imageFrameReg,stackUsed,"Registered Image ",lng=s.nz))
        stackUsed=s.stack_registered
        forwardButtPress(imageFrameReg,stackUsed,"Registered Image ")
 
########### Where most the calculations occur ##########
 
def calc():
    """
    Does the calculations that will happen every time 
    any parameter is changed or data is inputted
    """
    global threshold, maxpaths
    s.get_outliers(threshold=threshold,maxpaths=maxpaths)  
    s.get_averaged_image()
    global rijMaskFalse
    rijMaskFalse = s.show_Rij(mask=False,normalization=False,returnfig=True)
    global rijMaskTrue
    rijMaskTrue = s.show_Rij(mask=True,normalization=True,returnfig=True)
    global averageImage,fourierMask
    fourierMask=s.show_Fourier_mask_simple(returnfig=True)
    averageImage = s.show(returnfig=True)
    showView(root)
         
def loadData(fileType):
    """
    Imports data into program to be analyzed
    Also begins the analysis
    """
    t0=time()
    global s
    if fileType=='tif':
        try:
            file=askopenfilename()
            if file=='':
                return
            filepath=os.path.dirname(file)
            files=[]
            if(fileType=='tif'):
                for filename in listdir(filepath):
                    if not filename.startswith("."):
                        files.append(filepath+"/"+filename)
            im=np.array(Image.open(files[0]))
            fov=im.shape[0]
            stack=np.empty((fov,fov,len(files)))
            for i in range(len(files)):
                im=np.array(Image.open(files[i]))
                im=im/float(2**16)
                stack[:,:,i]=im	
        except:
            tkMessageBox.showwarning("Data Type Error","Please make sure you have selected the correct folder containing only .tif files")
            return
    elif fileType=='ser':
        
        try:
            file=askopenfilename()
            if file=='':
                return
            if not file.lower().endswith('.ser'):
                tkMessageBox.showwarning('Data Type Error',"Please make sure you have selected a .ser file to use this option.")
                return
            stack=serReader.serReader(file)	 
        except:
            tkMessageBox.showwarning('Data Type Error',"Please make sure you have selected a .ser file to use this option.")
            return
    else:
        try:
            file=askopenfilename()
            if file=='':
                return
            stack=np.rollaxis(imread(file),0,3)
            stack=stack[:,:,:]/float(2**16)
        except:
            tkMessageBox.showwarning('Data Type Error','Please make sure you have selected a .tif stack to use this option.')
            return
        
    # Instantiate imstack object
    s =stackregistration.imstack(stack)	 
    global findMaxima,loaded,fouriern,fourierMaskType
    s.getFFTs()
    s.makeFourierMask(mask=fourierMaskType,n=fouriern)
    s.findImageShifts(correlationType="cc",findMaxima=findMaxima,verbose=False)
    t=time()-t0
    s.set_nz(0,s.nz)
    loaded=True
    calc()
    
    print("Completed calculations in {} minutes {} seconds".format(int(t/60),t%60))
 
############ Redo Methods - Recalculate one part if parameters are changed ##########
   
def outRedo(bool,coef,z,shift):
    """
    Sets variables and redoes necessary calculations
    if outlier method parameters are changed
    """
    global polycoeffmax,zscore,shiftMax,outlierMethod
    polycoeffmax = int(coef)
    zscore = float(z)
    shiftMax = float(shift)
    if not bool:
        outlierMethod = "NN"
    else:
        outlierMethod="PF"
    global outpop
    outpop.destroy()
    calc()
    
def nzRedo(minText,maxText,listEntry):
    """
    Saves changes to the nz range and then redoes the calculations with new nz
    """
    s.set_nz(int(minText),int(maxText)+1)
    if listEntry!="":
        s.set_bad_images(map(int,listEntry.split(",")))
    global nzpop
    nzpop.destroy()
    calc()

def fourierRedo(n,type):
    """
    Updates fourier mask parameters and recalculates the mask
    Performs all other recalculations
    """
    global fouriern, fourierpop,fourierMaskType
    fouriern=int(n)
    fourierMaskType=type
    s.makeFourierMask(mask=type,n=fouriern)
    s.findImageShifts(correlationType="cc",findMaxima=findMaxima,verbose=False)
    calc()
    fourierpop.destroy()
 
def shiftRedo(maxima,numPeak,sigma,windowRad):
    """
    Updates changed parameters for finding the shift 
    Redoes all necessary calculations
    """
    global shiftpop,findMaxima, gaussiannumpeaks,sigmaguess,windowradius
    findMaxima=maxima
    gaussiannumpeaks=int(numPeak)
    sigmaguess=int(sigma)
    windowradius=int(windowRad)
    s.setGaussianFitParams(num_peaks=int(numPeak),sigma_guess=int(sigma),window_radius=int(windowRad))
    s.findImageShifts(correlationType="cc",findMaxima=maxima,verbose=False)
    shiftpop.destroy()
    calc()
          
########### Methods for Creating various pop-ups ##########
          
def outlierPopup():
    """
    Creates pop-up to edit the Outlier method
    """
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
        
        fitframe = tk.LabelFrame(outframe,text='Polynomial Fit Settings',relief='sunken',font=MEDIUMFONT)
        aroundframe = tk.LabelFrame(outframe,text='Nearest Points Fit Settings',relief='sunken',font=MEDIUMFONT)
        fitframe.grid(column=0,row=1,columnspan=2,sticky='W')
        aroundframe.grid(column=0,row=3,columnspan=2,sticky='W',pady=4)
        
        global polycoeffmax,zscore
        outCoeLabel=tk.Label(fitframe, text="Order of Polynomial fit:",font=SMALLFONT)
        outCoeLabel.grid(column=0,row=0)
        outCoeEntry = tk.Entry(fitframe, width=7, font=SMALLFONT)
        outCoeEntry.grid(column=1, row=0)
        outCoeEntry.insert(0,polycoeffmax)
        outZLabel = tk.Label(fitframe, text="Z Value:",font=SMALLFONT)
        outZLabel.grid(column=0, row=1)
        outZEntry = tk.Entry(fitframe, width=7, font=SMALLFONT)
        outZEntry.grid(column=1, row=1)
        outZEntry.insert(0,zscore)
        
        global shiftMax
        outShiftLabel=tk.Label(aroundframe,text="Max Shift:",font=SMALLFONT)
        outShiftLabel.grid(column=0,row=0)
        outShiftEntry = tk.Entry(aroundframe,width=7,font=SMALLFONT)
        outShiftEntry.grid(column=1,row=0)
        outShiftEntry.insert(0,shiftMax)
        
        outFitRadio=tk.Radiobutton(outframe, font=MEDIUMFONT,text="Use Polynomial Fit", variable=v, value=True,command=lambda: OutlierRadioUpdate(True,outCoeEntry,outZEntry,outShiftEntry))
        outAroundRadio=tk.Radiobutton(outframe, font=MEDIUMFONT, text="Use nearest points", variable=v, value=False,command=lambda: OutlierRadioUpdate(False,outCoeEntry,outZEntry,outShiftEntry))
        OutlierRadioUpdate(v.get(),outCoeEntry,outZEntry,outShiftEntry)
        outFitRadio.grid(column=0, row=0)
        outAroundRadio.grid(column=1,row=0)
        if v.get():
            outFitRadio.select()
        else:
            outAroundRadio.select()
        outSaveButton = tk.Button(outframe, text="Save", font=MEDIUMFONT,command=lambda: outRedo(v.get(),outCoeEntry.get(),outZEntry.get(),outShiftEntry.get()))
        outSaveButton.grid(column=0,row=4)
        outCancelButton = tk.Button(outframe, text="Cancel", font=MEDIUMFONT, command=outpop.destroy)
        outCancelButton.grid(column=1,row=4)
        
        outpop.mainloop()
    else:
       tkMessageBox.showwarning("Load Data","Please load images in to analyze first.") 
       
def nzPopup():
    """
    Creates pop-up to change the rangne of stack used
    """
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
        nzminLabel=tk.Label(nzframe, text="nz min:",font=SMALLFONT)
        nzminLabel.grid(column=0,row=0)
        nzminEntry = tk.Entry(nzframe, width=7, font=SMALLFONT)
        nzminEntry.grid(column=1, row=0)
        nzminEntry.insert(0,str(s.nz_min))
        nzmaxLabel = tk.Label(nzframe, text="nz max:",font=SMALLFONT)
        nzmaxLabel.grid(column=0, row=1)
        nzmaxEntry = tk.Entry(nzframe, width=7, font=SMALLFONT)
        nzmaxEntry.grid(column=1, row=1)
        nzmaxEntry.insert(0,str(s.nz_max-1))
        nzListLabel = tk.Label(nzframe,text="Bad images:",font=SMALLFONT)
        nzListLabel.grid(column=0,row=2)
        nzListEntry = tk.Entry(nzframe,width=7,font=SMALLFONT)
        nzListEntry.grid(column=1,row=2)
        if hasattr(s, 'bad_images'):
            nzListEntry.insert(0,','.join(map(str, s.bad_images)))
        nzSaveButton = tk.Button(nzframe, text="Save", font=MEDIUMFONT,command=lambda: nzRedo(nzminEntry.get(),nzmaxEntry.get(),nzListEntry.get()))
        nzSaveButton.grid(column=0,row=3)
        nzCancelButton = tk.Button(nzframe, text="Cancel", font=MEDIUMFONT, command=nzpop.destroy)
        nzCancelButton.grid(column=1,row=3)
    else:
        tkMessageBox.showwarning("Load Data","Please load images in to analze first")

def fourierPopup():
    """
    Creates pop-up to edit the Fourier Mask
    """
    global loaded
    if loaded:
        global fourierpop, fouriern,fourierMaskType
        fourierpop= tk.Tk()
        fourierpop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
        fourierpop.title("Edit Fourier Transform Mask")
        
        fourierExplainLabel=tk.Label(fourierpop,font=SMALLFONT,text="Top left: FFT of Slice 0\nwith mask shown\nTop right: FFT of slice 0\n with mask applied\nBottom: Cross\nCorrelation between\nslices 0 and "+str(int(s.nz/2)))
        fourierExplainLabel.grid(column=0,row=0,columnspan=2)
        fourierNLabel = tk.Label(fourierpop,text="N Value:",font=MEDIUMFONT)
        fourierNLabel.grid(column=0,row=1)

        fourierNEntry=tk.Entry(fourierpop,width=2,font=MEDIUMFONT)
        fourierNEntry.insert(0,fouriern)
        fourierNEntry.grid(column=1,row=1)
        
        radioText = tk.StringVar(fourierpop)
        radioText.set(fourierMaskType)
        fourierRadioFrame=tk.LabelFrame(fourierpop,text="Select Mask Type", relief="sunken",font=MEDIUMFONT)
        fourierNoneRadio=tk.Radiobutton(fourierRadioFrame,text="None",variable=radioText,value="none",font=SMALLFONT,command=lambda: fourierPopupDisplay("none",fourierNEntry.get(),fourierDisplayFrame))
        fourierLowRadio=tk.Radiobutton(fourierRadioFrame,text="Lowpass",variable=radioText,value="lowpass",font=SMALLFONT,command=lambda: fourierPopupDisplay("lowpass",fourierNEntry.get(),fourierDisplayFrame))
        fourierBandRadio=tk.Radiobutton(fourierRadioFrame,text="Bandpass",variable=radioText,value="bandpass",font=SMALLFONT,command=lambda: fourierPopupDisplay("bandpass",fourierNEntry.get(),fourierDisplayFrame))
        fourierHammingRadio=tk.Radiobutton(fourierRadioFrame,text="Hamming",variable=radioText,value="hamming",font=SMALLFONT,command=lambda: fourierPopupDisplay("hamming",fourierNEntry.get(),fourierDisplayFrame))
        fourierBlackmanRadio=tk.Radiobutton(fourierRadioFrame,text="Blackman",variable=radioText,value="blackman",font=SMALLFONT,command=lambda: fourierPopupDisplay("blackman",fourierNEntry.get(),fourierDisplayFrame))
        fourierGaussianRadio=tk.Radiobutton(fourierRadioFrame,text="Gaussian",variable=radioText,value="gaussian",font=SMALLFONT,command=lambda: fourierPopupDisplay("gaussian",fourierNEntry.get(),fourierDisplayFrame))
        fourierRadioFrame.grid(column=0,columnspan=2,row=2)
        fourierNoneRadio.grid(column=0,row=0)
        fourierLowRadio.grid(column=0,row=1)
        fourierBandRadio.grid(column=0,row=2)
        fourierHammingRadio.grid(column=0,row=3)
        fourierBlackmanRadio.grid(column=0,row=4)
        if radioText.get()=="none":
            fourierNoneRadio.select()
        elif radioText.get()=="lowpass" or radioText.get()=="hann":
            fourierLowRadio.select()
        elif radioText.get()=="bandpass":
            fourierBandRadio.select()
        elif radioText.get()=="hamming":
            fourierHammingRadio.select()
        elif radioText.get() == "blackman":
            fourierBlackmanRadio.select()
        fourierDisplayFrame=tk.Frame(fourierpop)
        fourierDisplayFrame.grid(column=3,row=0,rowspan=4)
        fourierPopupDisplay(fourierMaskType,fourierNEntry.get(),fourierDisplayFrame)
        fourierNEntry.bind("<Return>",lambda event, radio=radioText, n=fourierNEntry, frm=fourierDisplayFrame: fourierTraceUpdate(event,radio,n,frm))
        fourierNEntry.bind("<FocusOut>",lambda event, radio=radioText, n=fourierNEntry, frm=fourierDisplayFrame: fourierTraceUpdate(event,radio,n,frm))
        fourierButtonsFrame=tk.Frame(fourierpop)
        fourierButtonsFrame.grid(column=0,row=3)
        fourierSaveButton=tk.Button(fourierButtonsFrame,text="Save",font=MEDIUMFONT,command=lambda: fourierRedo(fourierNEntry.get(),radioText.get()))
        fourierSaveButton.grid(column=0,row=0)
        fourierCancelButton=tk.Button(fourierButtonsFrame,text="Cancel",font=MEDIUMFONT,command=fourierpop.destroy)
        fourierCancelButton.grid(column=1,row=0)
        fourierpop.mainloop()
    else:
        tkMessageBox.showwarning("Load Data","Please load images in to analyze first.")

def shiftPopup():
    """
    Creates pop-up to edit the method of finding shifts between images
    """
    global loaded
    if loaded:
        global shiftpop,findMaxima,gaussiannumpeaks,sigmaguess,windowradius,numiter,minwindowfrac
        shiftpop=tk.Tk()
        shiftpop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
        shiftpop.title("Edit Image Shift Calculation")    
        maximaFrame=tk.LabelFrame(shiftpop,text="Select Method for Finding Maxima",relief="ridge",font=MEDIUMFONT)
        maximaFrame.grid(row=0,column=0,columnspan=2)
        maximaRadio=tk.StringVar(shiftpop)
        maximaRadio.set(findMaxima)
        maximaPixelRadio=tk.Radiobutton(maximaFrame,text="Pixel",variable=maximaRadio,value="pixel",font=SMALLFONT,command=lambda: shiftRadioUpdate("pixel",numPeakEntry,sigmaEntry,windowRadEntry))
        maximaGfRadio=tk.Radiobutton(maximaFrame,text="Subpixel by Gaussian Fit",variable=maximaRadio,value="gf",font=SMALLFONT,command=lambda: shiftRadioUpdate("gf",numPeakEntry,sigmaEntry,windowRadEntry))
        maximaPixelRadio.grid(column=0,row=0)
        maximaGfRadio.grid(column=1,row=0)
        if findMaxima=="pixel":
            maximaPixelRadio.select()
        elif findMaxima=="gf":
            maximaGfRadio.select()
        else:
            maximaComRadio.select()    

        gfSettingFrame=tk.LabelFrame(shiftpop,text="Gaussian Fit Settings",relief="sunken",font=MEDIUMFONT)
        gfSettingFrame.grid(column=0,row=1,columnspan=2)
        numPeakLabel=tk.Label(gfSettingFrame,text="Number of Peaks:",font=SMALLFONT)
        numPeakEntry=tk.Entry(gfSettingFrame,width=2,font=SMALLFONT)
        numPeakEntry.insert(0,str(gaussiannumpeaks))
        sigmaLabel=tk.Label(gfSettingFrame,text="Sigma Guess:",font=SMALLFONT)
        sigmaEntry=tk.Entry(gfSettingFrame,width=2,font=SMALLFONT)
        sigmaEntry.insert(0,str(sigmaguess))
        windowRadLabel=tk.Label(gfSettingFrame,text="Window Radius:",font=SMALLFONT)
        windowRadEntry=tk.Entry(gfSettingFrame,width=2,font=SMALLFONT)
        windowRadEntry.insert(0,str(windowradius))
        numPeakLabel.grid(column=0,row=0)
        numPeakEntry.grid(column=1,row=0)
        sigmaLabel.grid(column=0,row=1)
        sigmaEntry.grid(column=1,row=1)
        windowRadLabel.grid(column=0,row=2)
        windowRadEntry.grid(column=1,row=2)
        shiftSaveButton=tk.Button(shiftpop,text="Save",font=MEDIUMFONT,command=lambda: shiftRedo(maximaRadio.get(),numPeakEntry.get(),sigmaEntry.get(),windowRadEntry.get(),))
        shiftSaveButton.grid(column=0,row=2)
        shiftCancelButton=tk.Button(shiftpop,text="Cancel",font=MEDIUMFONT,command=shiftpop.destroy)
        shiftCancelButton.grid(column=1,row=2)
        shiftRadioUpdate(maximaRadio.get(),numPeakEntry,sigmaEntry,windowRadEntry)
        shiftpop.mainloop()
    else:
        tkMessageBox.showwarning("Load Data","Please load images in to analyze first.") 
    
def aboutPopup():
    """
    Creates pop-up with information on this program
    """
    aboutPop=tk.Tk()
    aboutPop.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
    aboutPop.title("About")
    aboutPop.resizable(width=False, height=False)
    t1=tk.Label(aboutPop,text="Grr (Great Rigid Registration)",font=MEDIUMFONT)
    t1.grid(column=0,row=0)
    t2=tk.Label(aboutPop,text="v1.0   Released # Month 2017",font=SMALLFONT)
    t2.grid(column=0,row=1)
    t3=tk.Label(aboutPop,text="",font=SMALLFONT)
    t3.grid(column=0,row=2)
    t4=tk.Label(aboutPop,text=" Benjamin Savitzky, Emily Waite, Ismail El Baggari, Lena F. Kourkoutis",font=SMALLFONT)
    t4.grid(column=0,row=3)
    t5=tk.Label(aboutPop,text="Department of Applied and Engineering Physics, Cornell University",font=SMALLFONT)
    t5.grid(column=0,row=4)
 
########### Methods that help the pop-ups ###########
 
def fourierTraceUpdate(event, radio, n, frm):
    """
    Starts the figure update when the Enter key is pressed
    """
    fourierPopupDisplay(radio.get(), n.get(), frm) 

def fourierPopupDisplay(radio, n, frm):
    """
    Makes the figure in the fourier pop-up
    Updates every time something is changed
    """
    if n.isdigit():
        n=int(n)
        fft1=s.fftstack[:,:,0]
        fft2=s.fftstack[:,:,int(s.nz/2)]
        nx,ny = float(s.nx),float(s.ny)
        k_max=s.ny/n/2
        if radio=="bandpass":
            mask =np.fft.fftshift((s.kr<k_max)*(np.sin(2*n*np.pi*s.kr/ny)*np.sin(2*n*np.pi*s.kr/ny)))
        elif radio=="lowpass":
            mask = np.fft.fftshift((s.kr<k_max)*(np.cos(n*np.pi*s.kr/ny)**2))
        elif radio=="none":
            mask = np.ones_like(s.kr)
        elif radio=="blackman":
            mask = np.fft.fftshift((s.kr<k_max)*((21./50.)+0.5*np.cos((np.pi*s.kr)/k_max)+(2./25.)*np.cos((2*np.pi*s.kr)/k_max)))
        elif radio=="hamming":
            mask = np.fft.fftshift((s.kr<k_max)*((27./50.)+(23./50.)*np.cos((np.pi*s.kr)/k_max)))
        elif radio=="gaussian":
            mask= np.fft.fftshift(np.exp(-(s.kr/(k_max/3.))**2))
            
        cc = np.abs(np.fft.fftshift(s.fftw.ifft(mask * fft2 * np.conj(fft1))))
        gs=gridspec.GridSpec(3,2)
        gs.update(wspace=0.05,hspace=0.05)
        fig=plt.figure(figsize=(3,4.6))
        fig.subplots_adjust(wspace=0,hspace=0,left=0.01,right=0.99,top=0.99,bottom=0.01)
        ax1=plt.subplot(gs[0,0])
        ax2=plt.subplot(gs[0,1])
        ax=plt.subplot(gs[1:3,0:2])
        tempdetemp=np.log(np.abs(np.fft.fftshift(fft1)))
        ax1.matshow(tempdetemp[int(s.nx/4):int(3*s.nx/4),int(s.ny/4):int(3*s.ny/4)],cmap='gray',vmin=np.mean(tempdetemp.ravel()))
        ax1.matshow(np.fft.fftshift(mask)[int(s.nx/4):int(3*s.nx/4),int(s.ny/4):int(3*s.ny/4)],cmap='hot',alpha=0.3)
        ax2.matshow(np.log(np.abs(np.fft.fftshift(fft1*np.where(mask,mask,0.001))))[int(s.nx/4):int(3*s.nx/4),int(s.ny/4):int(3*s.ny/4)],cmap='gray',vmin=np.mean(np.log(np.abs(np.fft.fftshift(fft1)))[int(s.nx/4):int(3*s.nx/4),int(s.ny/4):int(3*s.ny/4)].ravel()))
        ax1.axis('off')
        ax2.axis('off')
        ax.matshow(cc,cmap="viridis",vmin=np.mean(np.log(np.abs(cc))).ravel())
        ax.axis('off')
        fourierCnvs=FigureCanvasTkAgg(fig, frm)
        fourierCnvs.show()
        fourierCnvs.get_tk_widget().grid(row=0,column=0)

def OutlierRadioUpdate(bool,outCoeEntry,outZEntry,outShiftEntry):
    """
    Makes the radio buttons in the outlier pop-up work properly
    """
    if bool:
        outCoeEntry.configure(state=tk.NORMAL)
        outZEntry.configure(state=tk.NORMAL)
        outShiftEntry.configure(state=tk.DISABLED)
    else:
        outCoeEntry.configure(state=tk.DISABLED)
        outZEntry.configure(state=tk.DISABLED)
        outShiftEntry.configure(state=tk.NORMAL)
   
def shiftRadioUpdate(selection,numPeakEntry,sigmaEntry,windowRadEntry):
    """
    Makes the radio buttons in the shift pop-up work properly
    """
    numPeakEntry.configure(state=tk.DISABLED)
    sigmaEntry.configure(state=tk.DISABLED)
    windowRadEntry.configure(state=tk.DISABLED)

    if selection=="gf":
        numPeakEntry.configure(state=tk.NORMAL)
        sigmaEntry.configure(state=tk.NORMAL)
        windowRadEntry.configure(state=tk.NORMAL)
                
########### Defines other menu items that are not pop-ups ############# 
 
def openHelp():
    """
    Opens User Guide pdf in user's default PDF program
    """
    filename=os.path.join(os.path.dirname(os.path.abspath(__file__)),"User_Guide.pdf")
    if sys.platform == "win32":
        os.startfile(filename)
    else:
        opener ="open" if sys.platform == "darwin" else "xdg-open"
        subprocess.call([opener, filename]) 

def saveReport():
    """
    Saves PDF report that includes averaged image, Rij matrices, FFT, etc. 
    """
    global loaded
    if loaded:
        file=tkFileDialog.asksaveasfilename(defaultextension='.pdf')
        if file=='':
            return
        s.save_report(file)
    else:
        tkMessageBox.showwarning("Load Data","Please load images in to analyze first.")

def saveImage():
    """
    Saves .tif of averaged image
    """
    global loaded
    if loaded:
        file=tkFileDialog.asksaveasfilename(defaultextension='.tif')
        if file=='':
            return
        s.save(file)
    else:
        tkMessageBox.showwarning("Load Data","Please load images in to analyze first.")

def saveTxt():
    """
    Saves the shift matrices to a .txt file
    """
    global loaded
    if loaded == False:
        tkMessageBox.showwarning("Load Data","Please load images in to analze first")
    else:
        fileName=tkFileDialog.asksaveasfilename(defaultextension='.txt')
        with file(fileName, 'w') as outfile:
            outfile.write('X shifts\n')
            np.savetxt(outfile, s.X_ij, fmt='%-7.2f')
            outfile.write('\nY Shifts\n')
            np.savetxt(outfile, s.Y_ij, fmt='%-7.2f')           
        
########## Root Set-up and Execution ##########
        
def rootSetup(root):
    """
    Makes the basic window with all of the menus
    Controls all of the menu options
    """
    root.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
    root.title('Grr')
    fig=plt.Figure()
    root.maxsize(width=int(11.5*fig.get_dpi()),height=int(7.5*fig.get_dpi()))
    root.resizable(width=False, height=False)
    menubar = tk.Menu(root)
    filemenu=tk.Menu(menubar, tearoff=0)
    filemenu.add_command(label="Import .tif stack",command=lambda:loadData(fileType='tifstack'))
    filemenu.add_command(label="Import .ser File",command=lambda: loadData(fileType='ser')) 
    filemenu.add_command(label="Import Extracted .tif files", command=lambda: loadData(fileType='tif'))
    filemenu.add_command(label="Save Image",command=saveImage)
    filemenu.add_command(label="Generate Report",command=saveReport )
    filemenu.add_command(label="Export Shift Data",command=saveTxt)
    editmenu=tk.Menu(menubar, tearoff=0)
    editmenu.add_command(label="Change nz range",command=nzPopup)
    editmenu.add_command(label="Change Outlier Method", command=outlierPopup)
    editmenu.add_command(label="Fourier Transformation Mask",command=fourierPopup)
    editmenu.add_command(label="Change Image Shift Method",command=shiftPopup) 
    helpmenu=tk.Menu(menubar, tearoff=0)
    helpmenu.add_command(label="About",command=aboutPopup)
    helpmenu.add_command(label="Getting Started",command=openHelp) 
    menubar.add_cascade(label="File", menu=filemenu)
    menubar.add_cascade(label="Edit", menu=editmenu)
    menubar.add_cascade(label="Help", menu=helpmenu)
    root.config(menu=menubar)        

loadParams()    
root = tk.Tk()
rootSetup(root)
root.protocol("WM_DELETE_WINDOW", quit)
root.mainloop()  
