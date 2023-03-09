from array import array
import os
import ctypes as c
import copy as cp
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from multiprocessing import shared_memory, Queue, Lock
from functools import partial
from scipy.ndimage import binary_dilation
from astropy.io import fits

def expansion(loop, erad):
    id = loop+1
    samp = white_arr_2d==id
    if erad>0:
        for k in range(erad):
            samp = binary_dilation(samp, structure = struc_arr_2d).astype(samp.dtype)
        eb = (white_arr_2d==0) * (samp!=0)
        white_arr_2d[eb] = cp.deepcopy(id)
    print('target :',id, '      ', end = '\r')

def to_numpy_array(white, shape_2d):
    white_arr_2d = np.ctypeslib.as_array(white)
    return white_arr_2d.reshape(shape_2d)

def init_worker(white, shape_2d, struc_shared):
    global white_arr_2d, struc_arr_2d
    white_arr_2d = to_numpy_array(white, shape_2d)
    struc_arr_2d = to_numpy_array(struc_shared, (3,3))

def mask_expansion(file_name):
    struc = np.zeros((3,3))
    xx = np.linspace(-1,1,3)
    nx, ny = np.meshgrid(xx,xx)
    nr = np.abs(nx) + np.abs(ny)
    struc[nr <= 1]= 1
    struc_shared = mp.Array(c.c_double, 9, lock = False)
    struc_shared[:] = np.ravel(struc)
    struc_arr_2d = to_numpy_array(struc_shared, (3,3))

    
    print(file_name)
    target, exten = file_name.split('.')
    segm = fits.open('%s_check.fits'%target)
    seg = segm[1].data

    os.system('mkdir expanded_map')

    if os.path.isfile('./%s_radius.txt'%target):
        rads = np.loadtxt('./%s_radius.txt'%target)
        print('    Load radius file')
    else:
        print('    No radius file. Make new one.')
        seg_max = np.int64(np.amax(seg))
        rads = np.zeros((seg_max))
        for i in range(seg_max):
            rads[i] = np.sqrt(np.sum(seg==i+1) / np.pi)
            print('        ',i+1, rads[i], end = '\r')
        print('    Made a new radius file for ', target)
        np.savetxt('%s_radius.txt'%target, rads)

    ces = np.linspace(0.1,5,50)

    x,y = seg.shape
    shape_2d = (x,y)
    white = mp.Array(c.c_double, x*y, lock = False)
    white[:] = np.ravel(seg)[:]
    white_arr_2d = to_numpy_array(white, shape_2d)

    for i in range(np.size(ces)):
        file_check = os.path.isfile('./expanded_map/%s_check_%1.1f.fits'%(target,ces[i]))
        if file_check == True:
            print('./expanded_map/%s_check_%1.1f.fits already exist. Skip!'%(target,ces[i]))
            continue
        try:
            segm_samp =  fits.open('./expanded_map/%s_check_%1.1f.fits'%(target,ces[i-1]))
            seg_samp = segm_samp[1].data
            seg_title = 'load previous ce'
            if i !=0:
                erads = np.int64(rads * ces[i]) - np.int64(rads*ces[i-1])
            else:
                erads = np.int64(rads * ces[i])
        except FileNotFoundError:
            seg_samp = cp.deepcopy(seg)
            erads = np.int64(rads * ces[i])
            seg_title = 'no saved previous ce'
        print('Exp cof : $.1f'%ces[i],end = '\r')

        queue = Queue()
        pool = mp.Pool(8, initializer=init_worker, initargs=(white,shape_2d,struc_arr_2d))
        pool.starmap(expansion, zip(np.arange(np.amax(seg)), erads))
        pool.close()
        pool.join()

        segm[1].data = white_arr_2d
        segm.writeto('./expanded_map/%s_check_%1.1f.fits'%(target,ces[i]), overwrite = True)

def sky_estimation(file_name):
    print(file_name, 'sky estimation start')
    target, exten = file_name.split('.')
    os.system('mkdir Orig_files')
    os.system('cp %s ./Orig_files/'%file_name)
    path2 = './expanded_map/'
    file_list2 = os.listdir(path2)
    mask_list = [file for file in file_list2 if file.startswith(target)]
    mask_list = np.sort(mask_list)
    mask_arr = np.linspace(0.1,0.1*np.size(mask_list),1*np.size(mask_list))
    dmask_arr = (mask_arr[:-1] + mask_arr[1:])/2

    if os.path.isfile('./skyregion.fits'):
        print('Load skyregion file.')
        skyregion = fits.open('skyregion.fits')
        skymap = skyregion[0].data
    else:
        print('No skyregion file. Make a new one.')
        skymap = np.zeros((2048,2048))
        for i in range(8):
            for j in range(8):
                skymap[2048//8*i:2048//8*(i+1),2048//8*j:2048//8*(j+1)] = j+i*8+1
        hdu = fits.PrimaryHDU(skymap)
        hdu.writeto('skyregion.fits', overwrite = True)

    sky_id = np.linspace(1,np.int64(np.amax(skymap)),np.int64(np.amax(skymap)))

    sci = fits.open(file_name)
    imag = sci[1].data
    j1 = 0
    j2 = 0
    sky_arr = np.zeros((sky_id.size,np.size(mask_list)))
    dsky_arr = np.zeros((sky_id.size,np.size(mask_list)-1))
    for f2 in mask_list:
        print('     ',f2, end = '\r')
        mask = fits.open(path2+f2)
        mask_map = mask[1].data
        for i in range(np.size(sky_id)):
            i2d_bool = (mask_map==0) * (skymap==i+1) * (imag!=0)
            sky_arr[i,j1] = np.nanmedian(imag[i2d_bool])
        j1 += 1
        try:
            dmask_map = mask_map - mask_map_old
            for i in range(np.size(sky_id)):
                i2d_bool = (dmask_map!=0) * (skymap==i+1) * (imag!=0)
                dsky_arr[i,j2] = np.nanmedian(imag[i2d_bool])
            j2 +=1
        except: NameError
        mask_map_old = cp.deepcopy(mask_map)
        mask.close()

    plt.figure(figsize = (6,4))
    ax1 = plt.subplot(111)
    ax2 = ax1.twinx()
    for i in range(np.size(sky_id)):
        ax1.plot(mask_arr, sky_arr[i,:], c = 'pink', alpha = 0.1)
    ax1.errorbar(mask_arr, np.nanmedian(sky_arr, axis = 0), yerr = np.nanstd(sky_arr, axis = 0), c = 'r', capsize = 0.5)
    ax2.plot(mask_arr, np.sum(np.isnan(sky_arr)==False, axis = 0), c = 'r', linestyle = '--')

    for i in range(np.size(sky_id)):
        ax1.plot(dmask_arr, dsky_arr[i,:], c = 'skyblue', alpha = 0.1)
    ax1.errorbar(dmask_arr, np.nanmedian(dsky_arr, axis = 0), yerr = np.nanstd(dsky_arr, axis = 0), c = 'b', capsize = 0.5)
    ax2.plot(dmask_arr, np.sum(np.isnan(dsky_arr)==False, axis = 0), c = 'b', linestyle = '--')
    ax1.set_xlabel(r'$c_{e}$')
    ax1.set_xlim(0,5.6)
    ax1.set_ylabel('SB')
    cond = np.sum(np.isnan(sky_arr)==False, axis = 0)>=50
    armin = np.argmin(np.nanmedian(sky_arr, axis = 0)[cond])
    ax1.axhline(np.nanmedian(sky_arr, axis = 0)[armin],c='k',linestyle = ':')
    ax1.axhline(np.nanmedian(sky_arr, axis = 0)[armin] + np.nanstd(sky_arr,axis = 0)[armin],c='k',linestyle = '-.')
    ax1.axhline(np.nanmedian(sky_arr, axis = 0)[armin] - np.nanstd(sky_arr,axis = 0)[armin],c='k',linestyle = '-.')
    ax1.set_ylim(np.nanmedian(sky_arr, axis = 0)[armin] - np.nanstd(sky_arr,axis = 0)[armin] * 10, np.nanmedian(sky_arr, axis = 0)[armin] + np.nanstd(sky_arr,axis = 0)[armin] * 10)
    ax2.set_ylabel('number of sky region')
    ax2.set_ylim(0,np.amax(skymap))
    plt.title(target)
    plt.savefig(target+'_sky.png', dpi = 200)
    plt.close()

    imag -= np.nanmedian(sky_arr, axis = 0)[armin]
    sci[1].data = imag
    sci.writeto(file_name, overwrite = True)

    os.system('rm '+path2+target+'*')


if __name__=='__main__':
    path = './'
    file_list = os.listdir(path)
    cal_list = [file for file in file_list if file.endswith("cal.fits")]
    for f in cal_list:
        mask_expansion(f)
        sky_estimation(f)

