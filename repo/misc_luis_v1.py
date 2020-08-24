# -*- coding: utf-8 -*-
""" 

Created on 26/06/19

Author : Carlos Eduardo Barbosa
Modified by Luis

"""
from __future__ import print_function, division

import os

import numpy as np
from astropy.table import Table
from astropy.io import fits

import context

def get_zps_dr1(Field, bands):
    """ Read the table containing the zero points for a given tile and given
    bands. """
    zpfile = os.path.join(context.tables_dir, "ZPfiles_2020_kadu",
                          "zps_tiles-FLUX_AUTO.fits")
    zpdata = fits.open(zpfile)
    zpdic = {a: {"R": b,
                 "F660": c,
                 "I": d}
              for a, b, c, d in zip(zpdata[1].data["TILE"],
                                    zpdata[1].data["ZP_R"],
                                    zpdata[1].data["ZP_F660"],
                                    zpdata[1].data["ZP_I"])} 
    zps = np.array([zpdic[Field][band] for band in bands])
    return zps
