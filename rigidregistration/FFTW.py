from __future__ import print_function
import sys, os
import numpy as np
if sys.version_info[0] < 3:
    import cPickle as pickle
else:
    import _pickle as pickle
relative_path = __file__
try:
    from pyfftw import *
    hasfftw = True
except ImportError:
    print('PYFFTW not found, please "pip install pyfftw" for up to 20x speedup')
    hasfftw = False


class WrapFFTW(object):
   def __init__(self, shape, **kwargs):
       self.shape = shape

       self._flags = kwargs.get('flags', ['FFTW_MEASURE'])
       self._threads = kwargs.get('threads', 8)
       self.wisdomfile = kwargs.get('wisdomfile', None)
       self.selfpath = os.path.split(relative_path)[0]

       if self.wisdomfile is None:
           self.wisdomfile = os.path.join(self.selfpath, 'fftw_wisdom.pkl')
       try:
           self._wisdom, self._gotwisdom = [], False
           with open(self.wisdomfile, 'rb') as infile:
               self._wisdom = pickle.load(infile)
           self._gotwisdom = import_wisdom(self._wisdom)
       except IOError:
           pass

       self.data = n_byte_align(np.zeros(self.shape), 16, 'complex128')
       self.data_k = n_byte_align(np.zeros(self.shape), 16, 'complex128')

       self.fft_object = FFTW(self.data, self.data_k,
                              axes=(0,1), flags = self._flags,
                              threads = self._threads)
       self.ifft_object = FFTW(self.data_k, self.data,
                               direction = 'FFTW_BACKWARD',
                               axes=(0,1), flags = self._flags,
                               threads = self._threads)

       self._wisdom = export_wisdom()
       with open(self.wisdomfile, 'wb') as outfile:
           pickle.dump(self._wisdom, outfile, -1)

   def fft(self, inp):
       self.data[:,:] = inp
       return self.fft_object().copy()

   def ifft(self, inp):
       self.data_k[:,:] = inp
       return self.ifft_object().copy()


class WrapFFTW_NUMPY(object):
   def __init__(self, shape, **kwargs):
       self.shape = shape

   def fft(self, inp):
       return np.fft.fftn(inp)

   def ifft(self, inp):
       return np.fft.ifftn(inp) 


if not hasfftw:
    WrapFFTW = WrapFFTW_NUMPY


