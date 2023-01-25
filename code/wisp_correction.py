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

def mad(a, mask = None, wisp_filt = None):
    delt = mask - a*wisp_filt
    delt = np.ravel(delt)
    delt = delt[np.isnan(delt)==False]
    mad = median_abs_deviation(delt)
    return(mad)

def wisp_correction(fname, output_extens = None):
    data = fits.open(fname)
    detector = data[0].header['DETECTOR'].lower()
    filter = data[0].header['FILTER']

    try: wisp_data = fits.open('/Volumes/Metal_Empire_001/JWST/wisp_templates/wisps_%s_%s.fits'%(detector,filter))
    except FileNotFoundError: wisp_data = None
    if wisp_data is None:
        print(fname, 'has no wisp features. Skip.')
        return None
    else:
        print(fname, 'has wisp features. Start the correction.')
        obs, order, point, chip, extens = fname.split('_')

        if (output_extens is None):
            try: os.mkdir('Orig_files')
            except FileExistsError: None
        old_fname = obs+'_'+order+'_'+point+'_'+chip+'_rate_before_wisp.fits'
        data.writeto('Orig_files/'+old_fname, overwrite = True)
        print(fname, 'flat start')

        ### Load files
        sci = data[1].data
        bool_0 = sci==0
        wisp = wisp_data[0].data

        ### Preprocess for wisp image
        wisp_filt = gaussian_filter(wisp, sigma = 2, mode = 'reflect')

        ### Preprocess for ramp image 
        bkg_estimator = MedianBackground()
        bkg = Background2D(sci, box_size = 64, filter_size = 7, bkg_estimator = bkg_estimator)
        threshold = 5.5 * bkg.background_rms

        segment_map = detect_sources(sci-bkg.background, threshold, npixels = 5)
        segmap = segment_map.data

        masked = cp.deepcopy(sci)
        masked[segmap!=0]=np.nan
        masked -= bkg.background_median

        ### Find the best a for MAD(maksed -  a * wisp_filt) to be minimized
        minimize_result = minimize(mad, 0, (masked, wisp_filt), method = 'Nelder-Mead', bounds = [[0,5]])

        ### Wisp correction
        sci -= minimize_result.x * wisp

        ### Saving
        sci[bool_0] = 0 
        data[1].data = sci
        if (output_extens is None):
            new_fname = cp.deepcopy(fname)
        else:
            new_fname = obs+'_'+order+'_'+point+'_'+chip+'_rate_'+output_extens+'.fits'
        data.writeto(new_fname, overwrite = True)
        data.close()
        print('saved as ', new_fname)

if __name__ == '__main__':
    path = './'
    file_list = os.listdir(path)

    rate_list = [file for file in file_list if file.endswith("rate.fits")]

    for f in rate_list:
        print(f)
        wisp_correction(f)