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
from scipy.ndimage import binary_dilation, median_filter

def tophat_filter(radii):
    struc = np.zeros((2*radii+1,2*radii+1))
    x = np.linspace(-radii,radii,2*radii+1)
    mx, my = np.meshgrid(x,x)
    mr = np.sqrt(mx**2 + my**2)
    struc[mr <= radii]= 1
    return struc

def ST(image, r_small):
    image_st_expan = median_filter(image, size = 7)
    struc = tophat_filter(r_small)
    image_small = binary_dilation(image_st_expan, structure = struc).astype(image_st_expan.dtype)
    return image_small

def LT(image, r_large):
    image_lt_expan = median_filter(image, size = 15)
    struc = tophat_filter(r_large)
    image_large = binary_dilation(image_lt_expan, structure = struc).astype(image_lt_expan.dtype)
    return image_large

def snowball_correction(fname, r_small=7, r_large=35, output_extens=None):
    # fname: input file name in str. the input extension should be '*ramp.fits'
    # r_small: radius for binary dilation of smaller tier (default = 7)
    # r_large: radius for binary dilation of larger tier (default = 35)
    # output_extens: output_extension for snowball corrected file. (default = None)
    #         if output_extens = None, it change the original file as 'Orig_files/*rate_old.fits' and save corrected image as '*rate.fits'
    data = fits.open(fname)
    obs, order, point, chip, exten = fname.split('_')

    if (output_extens is None):
        try: os.mkdir('Orig_files')
        except FileExistsError: None
        old_fname = obs+'_'+order+'_'+point+'_'+chip+'_ramp_old.fits'
        data.writeto('Orig_files/'+old_fname, overwrite = True)
    print(fname, 'snowball correction start.')

    ni, ng, nr, nc = np.shape(data[3].data)
    new_group_dq = np.zeros((ni, ng, nr, nc))

    for j in range(ni):
        for i in range(ng):
            group_dq = data[3].data[j,i,:,:]

            ### find where dq do not contain 4
            flags = bitmask.bitfield_to_boolean_mask(group_dq, ignore_flags = 4)
            flags[group_dq==0] = True

            ### delete other dq except 4
            group_dq_04 = cp.deepcopy(group_dq)
            group_dq_04[flags] = 0

            ### make smaller tier map and larger tier map and union them
            group_dq_st = ST(group_dq_04,r_small)
            group_dq_lt = LT(group_dq_04,r_large)
            group_dq_at = group_dq_st + group_dq_lt
            group_dq_at[group_dq_at!=0]=1
            group_dq_at[flags==False]=1

            ### from original dq, add united tier map
            group_dq[flags==False]-=4
            group_dq[group_dq_at!=0]+=4

            ### save new dq
            new_group_dq[j,i,:,:] = cp.deepcopy(group_dq)
            print(fname, ', %i of %i group_dqs are expanded.'%(j,i))
        
    data[3].data = cp.deepcopy(new_group_dq)
    if (output_extens is None):
        new_fname = cp.deepcopy(fname)
    else:
        new_fname = obs+'_'+order+'_'+point+'_'+chip+'_ramp_'+output_extens+'.fits'
    data.writeto(new_fname, overwrite = True)
    data.close()
    print('saved as ', new_fname)

if __name__=="__main__":
    path = './'
    file_list = os.listdir(path)
    uncal_list = [file for file in file_list if file.endswith("uncal.fits")]
    jyp = os.environ["JYP"]
    for f in uncal_list:
        os.system('strun '+jyp+'/code/detector1_uncal.params %s'%f)

    
    file_list = os.listdir(path)
    ramp_list = [file for file in file_list if file.endswith("ramp.fits")]

    for f in ramp_list:
        snowball_correction(f)
        os.system('strun '+jyp+'/code/detector1_ramp.params %s'%f)
