# -*- coding: utf-8 -*-
""" 
Created on 26/06/19

Author : Carlos Eduardo Barbosa
Note: Modified (improved) by Luis A. Gutiérrez 
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
import argparse
import sys

import context
import misc_luis_v1 as misc
from halpha_estimators import halpha_3F_method

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

parser = argparse.ArgumentParser(
    description="""Estimate weak Halpha emission if it exists""")

parser.add_argument("field", type=str,
                    default="STRIPE82-0164",
                    help="Name crop field centred in the source, taken the prefix (only the tile, ??_swp-crop.fits)")
parser.add_argument("--source", type=str,
                    default="39866",
                    help="Number source")

if __name__ == "__main__":
    cmd_args = parser.parse_args()
    tile = cmd_args.field
    Obj = cmd_args.source
    img_dir = os.path.join(context.data_dir, "Fields_Luis/"+tile)
    #tile = "STRIPE82-0164"
    bands = ['R', 'F660', 'I']
    filenames = [os.path.join(img_dir, "{}_{}_{}_swp-crop.fits".format(tile,
                                                                       band, Obj)) for band in bands]

    #hdu = fits.open(filenames[0])[0]
    try:
        data = np.array([fits.getdata(_) for _ in filenames])
    except FileNotFoundError:
        filenames = [os.path.join(img_dir, "{}_{}_swp".format(tile,
                           band)) for band in bands]
        data = np.array([fits.getdata(_) for _ in filenames])
        
    zps = misc.get_zps_dr1(tile, bands) # Get zero points from DR1 tables
    
    data *= np.power(10, -0.4 * zps )[:, None, None] # Apply zero point
    idxs = [bands.index(band) for band in ["F660", "R", "I"]]
    halpha = halpha_3F_method(data[idxs[0]], data[idxs[1]], data[idxs[2]]).value
    halpha_clip = np.clip(halpha, 0, np.infty)

    #halpha = np.clip(halpha, 0, np.infty)
    #plt.imshow(gaussian_filter(halpha, 0.2), origin="bottom")
    #plt.imshow(gaussian_filter(halpha, sigma=0.2), origin="bottom")
    #wcs = WCS(hdu.header)

    ####################################################################
    #Position##########################################################
    #####################################################################
    position = "position.reg"
    ra, dec = [], []

    try:
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
        
    except FileNotFoundError:
        print("File", position,
                "not found - is not necesary now")
    
    #PLOT
    #ax = plt.subplot(projection=wcs)#, label='overlays')
    f = plt.figure(figsize=(18,9))
    img = aplpy.FITSFigure(filenames[0], figure=f, subplot=(1, 2, 1), north=True)
    plt.imshow(gaussian_filter(halpha, 4.0), vmin=2.e3, vmax=3.e3, interpolation = 'nearest',
                                                                             origin='lower', cmap='gray', norm=LogNorm(), alpha=0.85)
    img.axis_labels.set_xtext('RA (J2000)')
    img.axis_labels.set_ytext('Dec (J2000)')
    img.axis_labels.set_font(size=20, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    img.tick_labels.set_font(size=20, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    
    img1 = aplpy.FITSFigure(filenames[0], figure=f, subplot=(1, 2, 2), north=True)
    plt.imshow(gaussian_filter(halpha, 2.0),  vmin=2.e3, vmax=3.e3, interpolation = 'nearest',
                                                          origin='lower', cmap=plt.cm.get_cmap('gray'), norm=LogNorm(), alpha=0.85)
    img1.axis_labels.set_xtext('RA (J2000)')
    img1.axis_labels.set_ytext('Dec (J2000)')
    img1.axis_labels.set_font(size=20, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')
    img1.tick_labels.set_font(size=20, weight='medium', stretch='normal', family='sans-serif', style='normal', variant='normal')

    img1.add_scalebar(10.0/60.)
    img1.scalebar.set_label('10 arcmin')
    img1.scalebar.set(color='blue', linewidth=8, alpha=1)
    img1.scalebar.set_font(size=30, weight='bold',
                      stretch='normal', family='sans-serif',
                      style='normal', variant='normal')
    img1.scalebar.set_linestyle('solid')
    #plt.colorbar()
   

    #img1.show_regions(os.path.join(img_dir,'position.reg'))
    #img1.show_markers(ra, dec, layer='marker', edgecolor='red', facecolor='none',  marker="o", s=10, alpha=0.35, linewidths=310.)#, layer='marker_set_1', edgecolor='black', facecolor='none', s=30, alpha=0.5, linewidths=20)

    img1.add_label(0.1, 0.93, "teste", color="white",
              horizontalalignment='left',
              weight='bold', size=25, relative=True, zorder=1000)
    dx, dy = 0.001, -0.001
    img1.add_label(0.1+dx, 0.93+dy, "teste", color="black", alpha=0.6,
              horizontalalignment='left',
              bbox={"facecolor": "black", "edgecolor": "none",# "pad": 20,
                    "alpha": 0.6, "boxstyle": "round,pad=0.5"},
              weight='bold', size=25, relative=True, zorder=999)

    img1.axis_labels.hide_y()
    img1.tick_labels.hide_y()
    
    # overlay = ax.get_coords_overlay('fk5')
    # overlay[0].set_axislabel('Right Ascension (J2000)')
    # overlay[1].set_axislabel('Declination (J2000)')
    # plt.colorbar()
    #ax.set_xlabel('Right Ascension')
    #ax.set_ylabel('Declination')
    img1.set_theme('publication')
    plt.savefig("halpha-{}_{}-ZP_kadu.pdf".format(cmd_args.field, 
                                          cmd_args.source))
