#!/usr/bin/env python

# Import global libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from time import time
from PIL import Image
from os import listdir
import subprocess
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
shiftText=None
averageImage=None
outpop = None
nzpop = None
loaded=False
outCoeEntry = None
outZEntry = None
outShiftEntry = None
slice=-1
v=None
LARGE_FONT = ("Verdana",24)
MEDIUMLARGE = ("Verdana",20)

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
            self.maskData=np.copy(s.mask)
            self.original=np.copy(s.mask)
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
        s.mask=self.maskData
        self.pop.destroy()
        recalc()
    def cancel(self):
        self.pop.destroy()
    def OnMouseDown(self,event):
        if sys.platform == "win32":
            xmin1=82
            xmax1=308
            ymin=63
            ymax=289
            xmin2=357
            xmax2=584
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
                s.mask=np.copy(self.maskData)
                self.img=s.show_Rij(mask=True,normalization=True)
                s.mask=np.copy(self.original)
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
                s.mask=np.copy(self.maskData)
                self.img=s.show_Rij(mask=True,normalization=True)
                s.mask=np.copy(self.original)
                self.drawCanvas()
    def drawCanvas(self):
        canvas = FigureCanvasTkAgg(self.img,self.frame)
        canvas.show()
        canvas.get_tk_widget().grid(column=0,row=0)
        canvas.get_tk_widget().config(highlightthickness=0)
        canvas.get_tk_widget().bind("<Button-1>", self.OnMouseDown)

#Displays the plots based on what is selected in View (Once view menu works, right now it just uses preselected if it works or not)
def showView(*args):
    currentx=0
    currenty=0
    global loaded
    if loaded:
        global averageImage
        figFrame=tk.Frame(root)
        figFrame.grid(column=0,row=0,pady=5,padx=5)
        canvas = FigureCanvasTkAgg(averageImage, figFrame)
        tools=tkagg.NavigationToolbar2TkAgg(canvas, figFrame)
        tools.update()
        canvas.show()
        canvas.get_tk_widget().config(highlightthickness=0)
        canvas.get_tk_widget().pack()
        tools.pack()

        outerframe=tk.Frame(root)
        imageFrame=tk.Frame(outerframe)
        backButton=tk.Button(outerframe, text="<", font=MEDIUMLARGE,command=lambda: backButtPress(imageFrame))
        forwardButton=tk.Button(outerframe,text=">",font=MEDIUMLARGE,command=lambda: forwardButtPress(imageFrame))
        backButton.pack(side="left",fill="y")
        forwardButton.pack(side="right",fill="y")
        imageFrame.pack(side="top",fill="both")
        outerframe.grid(row=0,column=1)
        forwardButtPress(imageFrame)

        global rijMaskFalse
        canvas2 = FigureCanvasTkAgg(rijMaskFalse, master=root)
        canvas2.get_tk_widget().config(highlightthickness=0)
        canvas2.show()
        canvas2.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
        canvas2.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)
        canvas2.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)
        
        global rijMaskTrue
        canvas3 = FigureCanvasTkAgg(rijMaskTrue, master=root)
        canvas3.get_tk_widget().config(highlightthickness=0)
        canvas3.show()
        canvas3.get_tk_widget().grid(column=1, row=1,padx=5,pady=5)
        canvas3.get_tk_widget().bind("<ButtonPress-1>",rijTrueOnMouseDown)
        canvas3.get_tk_widget().bind("<ButtonRelease-1>",rijTrueOnMouseRelease)

#For moving forward in the slice viewer        
def forwardButtPress(frame):
    global slice
    if slice==s.nz-1:
        slice=0
    else: 
        slice=slice+1
    fig,ax=plt.subplots(figsize=(3.75,3.75),dpi=100)
    ax.matshow(s.stack[:,:,slice],cmap='gray')
    ax.axis('off')
    ax.set_title("Stack Image "+str(slice))
    can=FigureCanvasTkAgg(fig,frame)
    can.get_tk_widget().config(highlightthickness=0)
    can.get_tk_widget().grid(column=0,row=0)
    
#For moving backward in the slice viewer
def backButtPress(frame):
    global slice
    if slice==0:
        slice=s.nz-1
    else:
        slice=slice-1
    fig,ax=plt.subplots(figsize=(3.75,3.75),dpi=100)
    ax.matshow(s.stack[:,:,slice],cmap='gray')
    ax.axis('off')
    ax.set_title("Stack Image "+str(slice))
    can=FigureCanvasTkAgg(fig,frame)
    can.get_tk_widget().config(highlightthickness=0)
    can.get_tk_widget().grid(column=0,row=0)

#Mask=False Rij Mouse Press
def rijFalseOnMouseDown(event):
    canvas2 = FigureCanvasTkAgg(s.show_Rij(mask=True),master=root)
    canvas2.get_tk_widget().config(highlightthickness=0)
    canvas2.show()
    canvas2.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
    canvas2.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)
    canvas2.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)

#Mask=False Rij Mouse Release
def rijFalseOnMouseRelease(event):
    global rijMaskFalse
    canvas2 = FigureCanvasTkAgg(rijMaskFalse,master=root)
    canvas2.get_tk_widget().config(highlightthickness=0)
    canvas2.show()
    canvas2.get_tk_widget().grid(column=0, row=1,pady=5,padx=5)
    canvas2.get_tk_widget().bind("<ButtonPress-1>",rijFalseOnMouseDown)
    canvas2.get_tk_widget().bind("<ButtonRelease-1>",rijFalseOnMouseRelease)

#Mask=True Rij Mouse Press
def rijTrueOnMouseDown(event):
    canvas3 = FigureCanvasTkAgg(s.show_Rij(mask=True), master=root)
    canvas3.get_tk_widget().config(highlightthickness=0)
    canvas3.show()
    canvas3.get_tk_widget().grid(column=1, row=1,padx=5,pady=5)
    canvas3.get_tk_widget().bind("<ButtonRelease-1>",rijTrueOnMouseRelease)
    canvas3.get_tk_widget().bind("<ButtonPress-1>",rijTrueOnMouseDown)

#Mask=True Rij Mouse Release
def rijTrueOnMouseRelease(event):
    global rijMaskTrue
    canvas3 = FigureCanvasTkAgg(rijMaskTrue, master=root)
    canvas3.get_tk_widget().config(highlightthickness=0)
    canvas3.show()
    canvas3.get_tk_widget().grid(column=1, row=1,padx=5,pady=5)
    canvas3.get_tk_widget().bind("<ButtonRelease-1>",rijTrueOnMouseRelease)
    canvas3.get_tk_widget().bind("<ButtonPress-1>",rijTrueOnMouseDown)
    
#Recalculates data without redoing entire thing if parameter is changed. Saves time. Good thing
def recalc(*args):
    try:
        s.get_averaged_image()
        global rijMaskFalse
        rijMaskFalse = s.show_Rij()
        global rijMaskTrue
        rijMaskTrue = s.show_Rij(mask=True,normalization=True)
        global averageImage
        averageImage = s.get_images()
        global shiftText
        shiftText=s.prnt
        showView()
    except ValueError:
        pass
        
#Runs when you select your data set
def loadData(fileType):
    try:
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
        else:
            stack=serReader.serReader(askopenfilename())

        # Instantiate imstack object
        s =stackregister.imstack(stack[:,:,:])	  
        s.getAllCrossCorrelations('subpixel')
        t=time()-t0
        s.set_nz(0,s.nz)
        global loaded
        loaded=True
        s.set_outliers()
        recalc()
        filemenu.entryconfig(2,command=lambda: s.save())
        filemenu.entryconfig(3, command=lambda: s.generate_report())
        filemenu.entryconfig(4,command=saveTxt)
        print("Completed calculations in {} minutes {} seconds".format(int(t/60),t%60))
    except ErrorValue:
        print(ErrorValue)
        
#Saves the shift matrixes to a text file
def saveTxt(*args):
    try:
        fileName=tkFileDialog.asksaveasfilename(defaultextension='.txt')
        with file(fileName, 'w') as outfile:
            outfile.write('#X shifts\n')
            np.savetxt(outfile, s.X_ij_full, fmt='%-7.2f')
            outfile.write('\n#Y Shifts\n')
            np.savetxt(outfile, s.Y_ij_full, fmt='%-7.2f')
    except:
        global loaded
        if loaded == False:
            tkMessageBox.showwarning("Load Data","Please load images in to analze first")
        else:
            pass
    
#This redoes the necessary parts if you change the outlier method    
def outRedo(bool,coef,z,shift):
    s.polycoeffmax = int(coef)
    s.z_score = float(z)
    s.method = bool
    s.max_shift = shift
    global outpop
    outpop.destroy()
    s.set_outliers()
    recalc()
    
#This redoes the necessary parts if you change the nz range
def nzRedo(minText,maxText):
    s.set_nz(int(minText),int(maxText))
    global nzpop
    nzpop.destroy()
    s.set_outliers()
    recalc()

#This makes the radio buttons in the change outlier method work
def radioUpdate(bool):
    global outCoeEntry
    global outZEntry
    global outShiftEntry
    global v
    v.set(bool)
    if bool:
        outCoeEntry.configure(state=tk.NORMAL)
        outZEntry.configure(state=tk.NORMAL)
        outShiftEntry.configure(state=tk.DISABLED)
    else:
        outCoeEntry.configure(state=tk.DISABLED)
        outZEntry.configure(state=tk.DISABLED)
        outShiftEntry.configure(state=tk.NORMAL)
        
#This creates the popup to change the settings when you change outlier method
def outlierPopup(*args):
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
        
        global v
        v.set(s.method)
        
        fitframe = tk.LabelFrame(outframe,text='Polynomial Fit Settings',relief='sunken')
        aroundframe = tk.LabelFrame(outframe,text='Nearest Points Fit Settings',relief='sunken')
        fitframe.grid(column=0,row=1,columnspan=2)
        aroundframe.grid(column=0,row=2,columnspan=2)
        
        global outCoeEntry
        global outZEntry
        outCoeLabel=tk.Label(fitframe, text="Coefficent of Polynomial fit:",font=LARGE_FONT)
        outCoeLabel.grid(column=0,row=0)
        outCoeEntry = tk.Entry(fitframe, width=7, font=LARGE_FONT)
        outCoeEntry.grid(column=1, row=0)
        outCoeEntry.insert(0,s.polycoeffmax)
        outZLabel = tk.Label(fitframe, text="Z Value:",font=LARGE_FONT)
        outZLabel.grid(column=0, row=1)
        outZEntry = tk.Entry(fitframe, width=7, font=LARGE_FONT)
        outZEntry.grid(column=1, row=1)
        outZEntry.insert(0,s.z_score)
        
        global outShiftEntry
        outShiftLabel=tk.Label(aroundframe,text="Max Shift:",font=LARGE_FONT)
        outShiftLabel.grid(column=0,row=0)
        outShiftEntry = tk.Entry(aroundframe,width=7,font=LARGE_FONT)
        outShiftEntry.grid(column=1,row=0)
        outShiftEntry.insert(0,s.max_shift)
        
        outFitRadio=tk.Radiobutton(outframe, text="Use Polynomial Fit", variable=v, value=True,command=lambda: radioUpdate(True))
        outAroundRadio=tk.Radiobutton(outframe, text="Use nearest points", variable=v, value=False,command=lambda: radioUpdate(False))
        radioUpdate(v.get())
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
       tkMessageBox.showwarning("Load Data","Please load images in to analze first") 
       
#This makes the popup to change the settings for the nz range
def nzPopup(*args):
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
        nzminEntry.insert(0,s.get_nzmin())
        nzmaxLabel = tk.Label(nzframe, text="nz max:",font=LARGE_FONT)
        nzmaxLabel.grid(column=0, row=1)
        nzmaxEntry = tk.Entry(nzframe, width=7, font=LARGE_FONT)
        nzmaxEntry.grid(column=1, row=1)
        nzmaxEntry.insert(0,s.get_nzmax())
        nzSaveButton = tk.Button(nzframe, text="Save", font=MEDIUMLARGE,command=lambda: nzRedo(nzminEntry.get(),nzmaxEntry.get()))
        nzSaveButton.grid(column=0,row=2)
        nzCancelButton = tk.Button(nzframe, text="Cancel", font=MEDIUMLARGE, command=nzpop.destroy)
        nzCancelButton.grid(column=1,row=2)
    else:
        tkMessageBox.showwarning("Load Data","Please load images in to analze first")

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
def openHelp():
    filename=os.path.join(os.path.dirname(os.path.abspath(__file__)),"User Guide.pdf")
    if sys.platform == "win32":
        os.startfile(filename)
    else:
        opener ="open" if sys.platform == "darwin" else "xdg-open"
        subprocess.call([opener, filename]) 
#Builds the menu bar and displays the window
root = tk.Tk()
root.iconbitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),"favicon.ico"))
root.title('Grr')
root.resizable(width=False, height=False)
menubar = tk.Menu(root)
v=tk.BooleanVar()

filemenu=tk.Menu(menubar, tearoff=0)
filemenu.add_command(label="Import .tif Files", command=lambda: loadData(fileType='tif'))
filemenu.add_command(label="Import .ser File",command=lambda: loadData(fileType='ser')) 
filemenu.add_command(label="Save Image")
filemenu.add_command(label="Generate Report")
filemenu.add_command(label="Export Shift Data")
editmenu=tk.Menu(menubar, tearoff=0)
editmenu.add_command(label="Change nz range",command=nzPopup)
editmenu.add_command(label="Change Outlier Method", command=outlierPopup)
editmenu.add_command(label="Custom mask",command=lambda: CustomMask(master=root)) 
helpmenu=tk.Menu(menubar, tearoff=0)
helpmenu.add_command(label="About",command=aboutPopup)
helpmenu.add_command(label="Getting Started",command=openHelp) 
menubar.add_cascade(label="File", menu=filemenu)
menubar.add_cascade(label="Edit", menu=editmenu)
menubar.add_cascade(label="Help", menu=helpmenu)

root.config(menu=menubar)
root.mainloop()  