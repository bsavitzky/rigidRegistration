# Image registration for tricky data

The rigidregistration package is designed for registering and averaging scanning transmission electron microscopy data, particulary in cases where low signal-to-noise ratios or periodicity-induced artifacts make registration difficult.


## A very quick overview

Aligning and averaging many images of a single sample region can vastly improve the quality of the final, averaged image.
However, in many applications, noisy data is unavoidal due to experimental considerations, and the resulting low SNRs may complicate or lead to incorrect image registration.
The basic ideas in this image registration package are as follows:

1. In a stack of images, the relative offsets between all pairs of images are calculated.
2. Physically, the relative offsets between different images must related to one another - e.g. the offsets between the first and second image and that between the second and third image should add up to the measured offset between the first and third image.
3. Combining the physical contraint imposed by (2) with the complete set of relative shifts from (1) allows a stack of images to be registered and averaged, even in the presence of (many) low-SNR induced errors.

For more detailed discussion, please see the publication associated with this package:
"Image registration of low signal-to-noise cryo-STEM data", Ultramicroscopy (2018), DOI: 10.1016/j.ultramic.2018.04.008.


## Getting started

The code can be acquired by cloning this repository to your computer, using the green "Clone or download" button, or by typing into the command line

```
git clone https://github.com/bsavitzky/rigidRegistration.git
```

Next, install the code by running the setup.py script. Navigate to the rigidregistration directory then type into the command line

```
python setup.py install
```

#### Dependencies

The dependencies of the package are: numpy, matplotlib, jupyter, and tifffile.
The package is built to optionally use pyfftw, which significantly speeds up calculations.
The first three dependencies are already available in most python implementations.
The last two can be installed from the command line with

```
pip install tifffile
pip install pyfftw
```


## Demo notebooks

The easiest way to familiarize yourself with the code, and get started using it on your own data, is with the demo notebooks found in the samplefiles directory.
These all .ipynb files, built to run from the jupyter notebook.
Launch jupyter from the command line with

```
jupyter notebook
```

There are three sample notebooks, described below.

1. DetailedWalkThrough.ipynb

This notebook is the best place to begin.
It works through registration of an example dataset, detailing exactly what is being done at each step of the process.
In doing so, it discusses all the basic features of the package that may be useful in handling your own data.

2. QuickWalkThrough.ipynb

This notebook works through registration of an example dataset, skipping over many of the details and focusing on obtaining a final, averaged image quickly.

3. MinimumWorkingExample.ipynb

This notebook provides the minimal code required to obtain a registered and averaged image.
Beginning here is not recommended, as very little explanation or comments are provided, including discussion of parameter selection which may be important for successful registration of your own data.
The code here may provide a useful starting point for automated batch processing.


## Thank you!

We hope you find this package useful in advancing your own research.
If this code is helpful to you and your work, please consider citing the associated publication:
"Image registration of low signal-to-noise cryo-STEM data", Ultramicroscopy (2018), DOI: 10.1016/j.ultramic.2018.04.008.




## Versioning

v. 1.0 - Original Release

## Authors

* **Benjamin H. Savitzky**
* **Emily N. Waite**
* **Lena F. Kourkoutis**

## License

This project is licensed under the MIT License.

## Acknowledgments

Many thanks to everyone who has been involved with this project at every level.
Particular thanks to the Kourkoutis research group at Cornell University for testing various versions of this code, finding bugs, and recommending essential improvements.
Thanks to Robert Hovden for the initial inspiration and assistance in getting this project started.
Thanks to Colin Clement for many useful discussions.
Thanks to everyone who provided samples and experimental data, including Ismail El Baggari, Berit H. Goodge, David J. Baek, John P. Sheckelton, Christopher Pasco, Hari Nair, Nathaniel J. Schreiber, Jason Hoffman, Alemayehu S. Admasu, Jaewook Kim, Sang-Wook Cheong, Anand Bhattacharya, Darrell G. Schlom, and Tyrel M. McQueen.
And thanks to you for using our code - we sincerely hope it is of use!


