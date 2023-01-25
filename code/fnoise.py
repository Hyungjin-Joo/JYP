import os
import copy as cp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clip
from photutils.segmentation import detect_sources
from photutils.background import Background2D, MedianBackground
from astropy.convolution import convolve, Gaussian2DKernel
from scipy.ndimage import binary_dilation


def fnoise_reduction(fname, bsize=None, fsize=None, radius=15, output_extens=None):
    # fname: input file name in str. the input extension should be '*rate.fits'
    # bsize: box size for Background2D (default = None)
    # fsize: filter size for Background2D (default = None)
    # radius: radius for binary dilation (default = 7)
    # output_extens: output_extension for 1/f noise reduced file. (default = None)
    #         if output_extens = None, it change the original file as 'Orig_files/*rate_old.fits' and save corrected image as '*rate.fits'


    data = fits.open(fname)
    obs, order, point, chip, exten = fname.split('_')

    if (output_extens is None):
        try: os.mkdir('Orig_files')
        except FileExistsError: None
        old_fname = obs+'_'+order+'_'+point+'_'+chip+'_rate_before_fnoise.fits'
        data.writeto('Orig_files/'+old_fname, overwrite = True)
    print(fname, '1/f noise reduction start.')

    tchip = ['nrca1','nrca2','nrca3','nrca4','nrcb1','nrcb2','nrcb3','nrcb4']

    if (bsize is None):
        if np.sum(np.isin(tchip, [chip]))!=1:
            print(fname, 'is LW. Set the boxsize = 32 for Background2D')
            bsize = 32
        else:
            print(fname, 'is SW. Set the boxsize = 64 for Background2D')
            bsize = 64

    if (fsize is None):
        if np.sum(np.isin(tchip, [chip]))!=1:
            print(fname, 'is LW. Set the filter size = 7 for Background2D')
            fsize = 7
        else:
            print(fname, 'is SW. Set the filter size = 13 for Background2D')
            fsize = 13

    sci = data[1].data
    cache = cp.deepcopy(sci)
    bool_0 = cache==0

    ### Background subtraction for saving ICL signals
    kernel = Gaussian2DKernel(x_stddev = 3.0, y_stddev = 3.0, x_size=5, y_size=5)  # FWHM = 3.
    convolved_sci = convolve(cache, kernel)

    bkg_estimator = MedianBackground()
    bkg = Background2D(convolved_sci, box_size = bsize, filter_size = fsize, bkg_estimator=bkg_estimator)

    sci -= bkg.background
    print(fname,', background subtraction done.')

    ### Segmentation map
    bkg2 = Background2D(sci, box_size = bsize, filter_size = fsize, bkg_estimator=bkg_estimator)
    threshold = 2.5 * bkg2.background_rms

    segment_map = detect_sources(sci, threshold, npixels=5)

    segmap = cp.deepcopy(segment_map.data)
    segmap[segmap !=0] = 1

    ### Segmap expansion with binary dilation
    struc1 = np.zeros((radius*2+1,radius*2+1))
    xx = np.linspace(-radius,radius,radius*2+1)
    nx, ny = np.meshgrid(xx,xx)
    nr = np.sqrt(nx**2 + ny**2)
    struc1[nr <= radius]= 1

    segmap = binary_dilation(segmap, structure=struc1).astype(segmap.dtype)
    print(fname,', mask map creation done.')

    ### Source masking
    # sci = cp.deepcopy(cache)
    sci[segmap!=0] = np.nan       

    ### 1/f noise estimation
    sci_hsplit = np.hsplit(sci, 4)
    row_noise = np.zeros((2048,2048))
    for i in range(4):
        samp_sci = sci_hsplit[i]
        filtered_sci = sigma_clip(samp_sci, sigma=2, maxiters=3, cenfunc = np.nanmedian, stdfunc = np.nanstd, axis = 1, masked = False)
        print(filtered_sci.shape)
        sci_hsplit_med = np.nanmedian(filtered_sci, axis = 1)
        sci_hsplit_med[np.isnan(sci_hsplit_med)]=0
        row_noise_hsplit = np.repeat(sci_hsplit_med.reshape(2048,1), 512, axis = 1)
        row_noise[0:2048,512*i:512*(i+1)] = cp.deepcopy(row_noise_hsplit)

    sci -= row_noise
    filtered_sci = sigma_clip(sci, sigma=2, maxiters=3, cenfunc = np.nanmedian, stdfunc = np.nanstd, axis = 0, masked = False)
    col_noise = np.nanmedian(filtered_sci, axis = 0)
    col_noise[np.isnan(col_noise)]=0
    col_noise = np.repeat(col_noise.reshape(1,2048), 2048, axis = 0)

    ### Noise reduction
    sci = cp.deepcopy(cache)
    sci -= row_noise
    sci -= col_noise
    print(fname,', 1/f noise correction done.')

    ### Saving
    sci[bool_0]=0
    data[1].data = sci
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
        fnoise_reduction(f)