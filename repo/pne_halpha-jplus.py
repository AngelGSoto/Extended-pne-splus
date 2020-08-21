# -*- coding: utf-8 -*-
""" 

Created on 26/06/19

Author : Carlos Eduardo Barbosa
Note: Modified (improved) by Luis A. Gutirérez 
      28/01/2020

"""
from __future__ import print_function, division

import os

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from astropy import coordinates as coord
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LogNorm
import aplpy

import context
import misc
from halpha_estimators import halpha_3F_method


if __name__ == "__main__":
    img_dir = os.path.join(context.data_dir, "jacoby")
    tile = "1000001-JPLUS-02363-v2"
    bands = ['rSDSS', 'J0660', 'iSDSS']
    filenames = [os.path.join(img_dir, "{}_{}_swp-65arc-crop.fits".format(tile,
                           band)) for band in bands]
    #hdu = fits.open(filenames[0])[0]

    data = np.array([fits.getdata(_) for _ in filenames])
    zps = misc.get_zps_dr1(tile, bands) # Get zero points from DR1 tables
    data *= np.power(10, -0.4 * zps )[:, None, None] # Apply zero point
    idxs = [bands.index(band) for band in ["J0660", "rSDSS", "iSDSS"]]
    halpha = halpha_3F_method(data[idxs[0]], data[idxs[1]], data[idxs[2]]).value
    print(data)
    #halpha = np.clip(halpha, 0, np.infty)
    #plt.imshow(gaussian_filter(halpha, 0.2), origin="bottom")
    #plt.imshow(gaussian_filter(halpha, sigma=0.2), origin="bottom")
    #wcs = WCS(hdu.header)

    ####################################################################
    #Position##########################################################
    #####################################################################
    position = "position.reg"
    ra, dec = [], []

    f = open(os.path.join(img_dir, position), 'r')
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()
        coor = line.split("(")[-1].split("\"")[0]
        ra1, dec1 = coor.split(",")[0:2]
        c = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg))
        ra.append(c.ra.degree)
        dec.append(c.dec.degree)
    
    #PLOT
    #ax = plt.subplot(projection=wcs)#, label='overlays')
    f = plt.figure(figsize=(18,9))
    img = aplpy.FITSFigure(filenames[0], figure=f, subplot=(1, 2, 1), north=True)
    plt.imshow(gaussian_filter(halpha, 3.),  vmin=2.e3, vmax=3.e3, origin="bottom", cmap="gray", norm=LogNorm(), alpha=0.75)
    img.axis_labels.set_xtext('RA (J2000)')
    img.axis_labels.set_ytext('Dec (J2000)')
    img.axis_labels.set_font(size=20, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    img.tick_labels.set_font(size=20, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    
    img1 = aplpy.FITSFigure(filenames[0], figure=f, subplot=(1, 2, 2), north=True)
    plt.imshow(gaussian_filter(halpha, 5.),  vmin=2.e3, vmax=3.e3, origin="bottom", cmap="gray", norm=LogNorm(), alpha=0.75)
    img1.axis_labels.set_xtext('RA (J2000)')
    img1.axis_labels.set_ytext('Dec (J2000)')
    img1.axis_labels.set_font(size=20, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    img1.tick_labels.set_font(size=20, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')

    img1.add_scalebar(8.0/60.)
    img1.scalebar.set_label('8 arcmin')
    img1.scalebar.set(color='blue', linewidth=8, alpha=1)
    img1.scalebar.set_font(size=30, weight='bold',
                      stretch='normal', family='sans-serif',
                      style='normal', variant='normal')
    img1.scalebar.set_linestyle('solid')
   

    #img1.show_regions(os.path.join(img_dir,'position.reg'))
    img1.show_markers(ra, dec, layer='marker', edgecolor='red', facecolor='none',  marker="o", s=10, alpha=0.35, linewidths=310.)#, layer='marker_set_1', edgecolor='black', facecolor='none', s=30, alpha=0.5, linewidths=20)

    img1.add_label(0.1, 0.93, "Jacoby 1", color="white",
              horizontalalignment='left',
              weight='bold', size=25, relative=True, zorder=1000)
    dx, dy = 0.001, -0.001
    img1.add_label(0.1+dx, 0.93+dy, "Jacoby 1", color="black", alpha=0.6,
              horizontalalignment='left',
              bbox={"facecolor": "black", "edgecolor": "none",# "pad": 20,
                    "alpha": 0.6, "boxstyle": "round,pad=0.5"},
              weight='bold', size=25, relative=True, zorder=999)

    img1.axis_labels.hide_y()
    img1.tick_labels.hide_y()
    
    # overlay = ax.get_coords_overlay('fk5')
    # overlay[0].set_axislabel('Right Ascension (J2000)')
    # overlay[1].set_axislabel('Declination (J2000)')
    #plt.colorbar()
    #ax.set_xlabel('Right Ascension')
    #ax.set_ylabel('Declination')
    img1.set_theme('publication')
    plt.savefig("halpha-jacoby-teste.pdf")
    #plt.show()

