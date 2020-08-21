# -*- coding: utf-8 -*-
""" 

Created on 26/06/19

Author : Carlos Eduardo Barbosa

"""
from __future__ import print_function, division

import os

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

import context
import misc
from halpha_estimators import halpha_3F_method
from matplotlib.colors import LogNorm


if __name__ == "__main__":
    img_dir = os.path.join(context.data_dir, "pne1")
    tile = "HYDRA-0054"
    #bands = ['U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660', 'I',
             #'F861', 'Z']
    bands = ['R', 'F660', 'I']
    filenames = [os.path.join(img_dir, "{}_{}_swp-crop.fits".format(tile,
                           band)) for band in bands]
    data = np.array([fits.getdata(_) for _ in filenames])
    zps = misc.get_zps_dr1(tile, bands) # Get zero points from DR1 tables
    data *= np.power(10, -0.4 * zps)[:, None, None] # Apply zero point
    idxs = [bands.index(band) for band in ["F660", "R", "I"]]
    halpha = halpha_3F_method(data[idxs[0]], data[idxs[1]], data[idxs[2]]).value
    halpha = np.clip(halpha, 0, np.infty)
    plt.imshow(gaussian_filter(halpha, 0.1), vmin=2.e3, vmax=3.e3, cmap="gray", origin="lower",  norm=LogNorm(), alpha=0.9)
    plt.colorbar()
    plt.show()

