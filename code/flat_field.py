import os
import copy as cp
from tokenize import group
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.nddata import bitmask
from astropy.stats import sigma_clip
from photutils.segmentation import detect_sources
from photutils.background import Background2D, MedianBackground
from astropy.convolution import convolve, Gaussian2DKernel
from scipy.ndimage import binary_dilation, median_filter, gaussian_filter
from scipy.stats import median_abs_deviation
from scipy.optimize import minimize

def flat_field(fname, output_extens = None):
    data = fits.open(fname)
    detector = data[0].header['DETECTOR'].lower()
    filter = data[0].header['FILTER']

    try: flat_data = fits.open('/Volumes/Metal_Empire_001/JWST/flat_templates/flat_%s_%s.fits'%(detector,filter))
    except FileNotFoundError: flat_data = None
    if flat_data is None:
        print(fname, 'has no flat data. Skip.')
        return None

    flat = flat_data[0].data

    obs, order, point, chip, exten = fname.split('_')

    if (output_extens is None):
        try: os.mkdir('Orig_files')
        except FileExistsError: None
        old_fname = obs+'_'+order+'_'+point+'_'+chip+'_rate_before_flat.fits'
        data.writeto('Orig_files/'+old_fname, overwrite = True)
    print(fname, 'flat start')

    imag = data[1].data
    imag /= flat

    imag[flat==0]=0

    data[1].data = imag
    if (output_extens is None):
        new_fname = cp.deepcopy(fname)
    else:
        new_fname = obs+'_'+order+'_'+point+'_'+chip+'_rate_'+output_extens+'.fits'
    data.writeto(new_fname, overwrite = True)
    data.close()
    print('saved as ', new_fname)

if __name__=="__main__":
    path = './'
    file_list = os.listdir(path)
    rate_list = [file for file in file_list if file.endswith("rate.fits")]

    for f in rate_list:
        flat_field(f)