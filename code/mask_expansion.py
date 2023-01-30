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

def mask_expansion():
    path = './'
    file_list = os.listdir(path)
    cal_list = [file for file in file_list if file.endswith("cal.fits")]

    struc = np.zeros((3,3))
    xx = np.linspace(-1,1,3)
    nx, ny = np.meshgrid(xx,xx)
    nr = np.abs(nx) + np.abs(ny)
    struc[nr <= 1]= 1
    struc_shared = mp.Array(c.c_double, 9, lock = False)
    struc_shared[:] = np.ravel(struc)
    struc_arr_2d = to_numpy_array(struc_shared, (3,3))

    for tid in cal_list:
        print(tid)
        target, exten = tid.split('.')
        segm = fits.open('%s_check.fits'%target)
        seg = segm[1].data

        os.system('mkdir expanded_map')

        if os.path.isfile('./%s_radius.txt'%target):
            rads = np.loadtxt('./%s_radius.txt'%target)
        else:
            seg_max = np.int64(np.amax(seg))
            rads = np.zeros((seg_max))
            for i in range(seg_max):
                rads[i] = np.sqrt(np.sum(seg==i+1) / np.pi)
                print(i+1, rads[i], end = '\r')
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
                print(ces[i],'skip!')
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
            print(seg_title)
            print('Exp cof :', ces[i], 'start')

            queue = Queue()
            pool = mp.Pool(8, initializer=init_worker, initargs=(white,shape_2d,struc_arr_2d))
            pool.starmap(expansion, zip(np.arange(np.amax(seg)), erads))
            pool.close()
            pool.join()

            segm[1].data = white_arr_2d
            segm.writeto('./expanded_map/%s_check_%1.1f.fits'%(target,ces[i]), overwrite = True)

