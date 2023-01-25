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

def tophat_filter(radii):
    struc = np.zeros((2*radii+1,2*radii+1))
    x = np.linspace(-2*radii+1,2*radii+1)
    mx, my = np.meshgrid(x,x)
    mr = np.sqrt(mx**2 + my**2)
    struc[mr <= radii]= 1
    return struc

def flat_maker(fnames, bsize=None, fsize=None, radius=7, output_extens=None):
    flats = np.zeros((np.size(fnames),2048,2048))
    i = 0
    path = './'
    file_list = os.listdir(path)

    mask_list = [file for file in file_list if file.endswith("mask.fits")]

    for fname in fnames:
        obs, order, point, chip, exten = fname.split('_')

        ### Check if there is maksed map
        mask_name = obs+'_'+order+'_'+point+'_'+chip+'_mask.fits'
        if np.sum(np.isin(mask_list, mask_name))!=0:
            print(fname, 'has masked image')
            data = fits.open(mask_name)
            sci = data[1].data
            flats[i,:,:]=sci
            i+=1
            data.close()
            continue

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
        

        data = fits.open(fname)
        sci = data[1].data
        cache = cp.deepcopy(sci)

        if (output_extens is None):
            try: os.mkdir('Orig_files')
            except FileExistsError: None
            old_fname = obs+'_'+order+'_'+point+'_'+chip+'_rate_old.fits'
            data.writeto('Orig_files/'+old_fname, overwrite = True)
        print(fname, i, '/', np.size(fnames))

        ### Background subtraction for saving ICL signals
        kernel = Gaussian2DKernel(x_stddev = 3.0, y_stddev = 3.0, x_size=5, y_size=5)  # FWHM = 3.
        convolved_sci = convolve(cache, kernel)

        bkg_estimator = MedianBackground()
        bkg = Background2D(convolved_sci, box_size = bsize, filter_size = fsize, bkg_estimator=bkg_estimator)

        sci -= bkg.background
        print(fname,', background subtraction done.')

        ### Segmentation map
        bkg2 = Background2D(sci, box_size = bsize, filter_size = fsize, bkg_estimator=bkg_estimator)
        threshold = 5.5 * bkg2.background_rms

        segment_map = detect_sources(sci, threshold, npixels=5)

        segmap = cp.deepcopy(segment_map.data)
        if  np.sum(np.isin(tchip, [chip])) != 1:
            struc1 = np.zeros((31,31))
            xx = np.linspace(-15,15,31)
            nx, ny = np.meshgrid(xx,xx)
            nr = np.sqrt(nx**2 + ny**2)
            struc1[nr <= 15]= 1
        else:
            struc1 = np.zeros((7,7))
            xx = np.linspace(-3,3,7)
            nx, ny = np.meshgrid(xx,xx)
            nr = np.sqrt(nx**2 + ny**2)
            struc1[nr <= 3]= 1
        segmap = binary_dilation(segmap, structure=struc1).astype(segmap.dtype)
        print(fname,', mask map creation done.')
        
        cache[segmap!=0] = np.nan
        data[1].data = cache
        data.writeto(obs+'_'+order+'_'+point+'_'+chip+'_'+'mask.fits', overwrite = True)
        flats[i,:,:] = cache/np.nanmedian(cache)
        data.close()
        i+=1
        print(fname,' done')
    
    # flat = sigma_clip(flats, sigma = 5, cenfunc = np.nanmedian, stdfunc = np.nanstd, masked = False, axis = 0)
    # plt.figure()
    # plt.imshow(flat[0,:,:])
    # plt.show()
    flat = np.nanmedian(flats, axis=0)
    flat[np.isnan(flat)]=0
    flat /= np.nanmedian(flat[flat!=0])
    flat /= np.nanmedian(flat[1022:1026,1022:1026])
    hdu = fits.PrimaryHDU(flat)
    hdu.writeto('flat.fits', overwrite = True)

def flat_applier(fname, flat, output_extens):
    data = fits.open(fname)
    obs, order, point, chip, exten = fname.split('_')

    if (output_extens is None):
        try: os.mkdir('Orig_files')
        except FileExistsError: None
        old_fname = obs+'_'+order+'_'+point+'_'+chip+'_rate_old.fits'
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

    flat_maker(rate_list)
    
    flat = fits.open('flat.fits')
    flat_imag = flat[0].data
    for f in rate_list:
        flat_applier(f,flat_imag, output_extens='flat')