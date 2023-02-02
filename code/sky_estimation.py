import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
import copy as cp
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import median_abs_deviation as mad

path = './'
file_list = os.listdir(path)
cal_list = [file for file in file_list if file.endswith("nrcb1_cal.fits")]

for tid in cal_list:
    print(tid)
    target, exten = tid.split('.')
    os.system('mkdir Orig_files')
    os.system('cp %s ./Orig_files/'%tid)
    path2 = './expanded_map/'
    file_list2 = os.listdir(path2)
    mask_list = [file for file in file_list2 if file.startswith(target)]
    mask_list = np.sort(mask_list)
    mask_arr = np.linspace(0.1,0.1*np.size(mask_list),1*np.size(mask_list))
    dmask_arr = (mask_arr[:-1] + mask_arr[1:])/2

    if os.path.isfile('./skyregion.fits'):
        skyregion = fits.open('skyregion.fits')
        skymap = skyregion[0].data
    else:
        skymap = np.zeros((2048,2048))
        for i in range(8):
            for j in range(8):
                skymap[2048//8*i:2048//8*(i+1),2048//8*j:2048//8*(j+1)] = j+i*8+1
        hdu = fits.PrimaryHDU(skymap)
        hdu.writeto('skyregion.fits', overwrite = True)

    sky_id = np.linspace(1,np.int64(np.amax(skymap)),np.int64(np.amax(skymap)))

    sci = fits.open(tid)
    imag = sci[1].data
    j1 = 0
    j2 = 0
    sky_arr = np.zeros((sky_id.size,np.size(mask_list)))
    dsky_arr = np.zeros((sky_id.size,np.size(mask_list)-1))
    for f2 in mask_list:
        print('     ',f2)
        mask = fits.open(path2+f2)
        mask_map = mask[1].data
        for i in range(np.size(sky_id)):
            i2d_bool = (mask_map==0) * (skymap==i+1) * (imag!=0)
            sky_arr[i,j1] = np.nanmedian(imag[i2d_bool])
            print('sky',i+1, np.nanmedian(imag[i2d_bool]), '                ', end = '\r')
        j1 += 1
        try:
            dmask_map = mask_map - mask_map_old
            for i in range(np.size(sky_id)):
                i2d_bool = (dmask_map!=0) * (skymap==i+1) * (imag!=0)
                dsky_arr[i,j2] = np.nanmedian(imag[i2d_bool])
                print('dsky',i+1, np.nanmedian(imag[i2d_bool]), '                ', end = '\r')
            j2 +=1
        except: NameError
        mask_map_old = cp.deepcopy(mask_map)
        mask.close()

    plt.figure(figsize = (6,4))
    ax1 = plt.subplot(111)
    ax2 = ax1.twinx()
    for i in range(np.size(sky_id)):
        ax1.plot(mask_arr, sky_arr[i,:], c = 'pink', alpha = 0.1)
    ax1.errorbar(mask_arr, np.nanmedian(sky_arr, axis = 0), yerr = mad(sky_arr[np.isnan(sky_arr)==False], axis = 0), c = 'r', capsize = 0.5)
    ax2.plot(mask_arr, np.sum(np.isnan(sky_arr)==False, axis = 0), c = 'r', linestyle = '--')

    for i in range(np.size(sky_id)):
        ax1.plot(dmask_arr, dsky_arr[i,:], c = 'skyblue', alpha = 0.1)
    ax1.errorbar(dmask_arr, np.nanmedian(dsky_arr, axis = 0), yerr = mad(dsky_arr[np.isnan(dsky_arr)==False], axis = 0), c = 'b', capsize = 0.5)
    ax2.plot(dmask_arr, np.sum(np.isnan(dsky_arr)==False, axis = 0), c = 'b', linestyle = '--')
    ax1.set_xlabel(r'$c_{e}$')
    ax1.set_xlim(0,5.6)
    ax1.set_ylabel('SB')
    cond = np.sum(np.isnan(sky_arr)==False, axis = 0)>=20
    armin = np.argmin(np.nanmedian(sky_arr, axis = 0)[cond])
    ax1.axhline(np.nanmedian(sky_arr, axis = 0)[armin],c='k',linestyle = ':')
    ax1.axhline(np.nanmedian(sky_arr, axis = 0)[armin] + mad(sky_arr[np.isnan(sky_arr)==False],axis = 0),c='k',linestyle = '-.')
    ax1.axhline(np.nanmedian(sky_arr, axis = 0)[armin] - mad(sky_arr[np.isnan(sky_arr)==False],axis = 0),c='k',linestyle = '-.')
    ax1.set_ylim(np.nanmedian(sky_arr, axis = 0)[armin] - mad(sky_arr[np.isnan(sky_arr)==False],axis = 0) * 10, np.nanmedian(sky_arr, axis = 0)[armin] + mad(sky_arr[np.isnan(sky_arr)==False],axis = 0) * 10)
    ax2.set_ylabel('number of sky region')
    ax2.set_ylim(0,np.amax(skymap))
    plt.title(target)
    plt.savefig(target+'_sky.png', dpi = 200)
    plt.close()

    imag -= np.nanmedian(sky_arr, axis = 0)[armin]
    sci[1].data = imag
    sci.writeto(tid, overwrite = True)
