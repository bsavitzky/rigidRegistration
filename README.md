# Rigid Registration

We have created a new and imporved method of registering stacks of images from scanning transmission electron microscopy (STEM).
In this method, we compare every pair of images, rather than comparing all images to one reference and use the redundant information from all of these pairs to determine the shifts between each image before averaging and getting the final image. 
See our paper in *not gonna jinx this* for more information on our methods. 

## Getting Started

To start, download the entire project and put all files in a folder somewhere in your file directory where you can point to later. 
Then make sure you have all of the prerequisites.

### Prerequisites

In order to run this program, you need either Python 2 or Python 3. If you do not already have python installed on your system, we recommend installing a version of [conda](https://conda.io/docs/user-guide/install/download.html).
Either Anaconda or Miniconda works well, but if you use Miniconda, you may have to install more packages. 
Once this is installed, check to make sure you have all the packages you need. The best way to do this is open your command prompt and run

```
pip install Jupyter
pip install Matplotlib
pip install Tifffile
pip install pillow
pip install pyfftw
```
Should you already have these packages, an error will appear saying that "Requirement already satisfied". Otherwise, the package will be installed. 

NOTE: pyfftw is optional. 


### Running

To use the python script, we have a few different implimations to use. 
If you are unfamiliar with the python language and/or registration, using the graphic user interphase may be best. 
To run, open a command prompt and enter


```
python C:\Path\To\Registration\Download\rigidregistration\gui.py
```

This will open a screen with a menu bar at the top. From the menu bar, go to File > Import .tif Stack.
From the file explorer opened, navigate to any data set stored in a .tif you have or navigate to our sample data set in rigidRegistration\samplefiles.
After opening this file, the data will be registered and the output will display most steps of the registation process. For more information on the graphic user interphase, see [the user guide](User Guide.pdf).

The other implimation to look at is in an iPython notebook. This is best if you are familiar with python already. To start this, run 

```
Jupter notebook
```

From the page that opens, navigate to the file directory where the rigid registration download is. In that folder go to samplefiles\sample_notebook.ipynb.
This code already points to the sample data set we provide. Run each cell in succession (Shift+Enter will run the current cell). Read the descriptions in the notebook as needed to understand what each cell does and how the different functions of our rigid registration method applies to varying data sets.


## Versioning

v. 1.0 - Original Release

## Authors

* **Benjamin H. Savitzky** - Original Developer
* **Emily Waite** - Graphic User Interphase Developer

## License

This project is licensed under the MIT License?

## Acknowledgments

* Stuff about how cool Cornell is
* And Professor Lena Kourkoutis :)
