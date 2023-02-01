import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import imshow_norm, ManualInterval, LogStretch

path = './'
file_list = os.listdir(path)
exten = input('Extension : ')
rate_list = [file for file in file_list if file.endswith(exten+".fits")]
try: os.mkdir('./figs')
except: FileExistsError

for f in rate_list:
    data = fits.open(f)
    sci = data[1].data
    print(f)
    plt.figure()
    # if np.nanmedian(sci)>0:
    #     plt.imshow(sci, origin = 'lower', vmin = -10* np.nanmedian(sci), vmax = 10.0*np.nanmedian(sci))
    # elif np.nanmedian(sci)==0:
    #     plt.imshow(sci, origin = 'lower', vmin = -0.05, vmax = 0.05)
    # else:
    #     plt.imshow(sci, origin = 'lower', vmin = 1.5* np.nanmedian(sci), vmax = -1.5*np.nanmedian(sci))
    # plt.imshow(sci, origin='lower', vmin = np.nanmedian(sci) - np.nanstd(sci), vmax = np.nanmedian(sci) + np.nanstd(sci))
    im, show = imshow_norm(sci, origin = 'lower', interval = ManualInterval(vmin = np.nanmedian(sci), vmax = np.nanmedian(sci) * 10), stretch = LogStretch(), cmap = 'Greys_r')
    plt.savefig('./figs/%s.png'%f)
    plt.close()