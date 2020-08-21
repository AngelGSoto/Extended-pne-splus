# -*- coding: utf-8 -*-
""" 

Created on 26/06/19

Author : Carlos Eduardo Barbosa

"""
from __future__ import print_function, division

import os

import numpy as np
from astropy.table import Table

import context

def get_zps_dr1(tile, bands):
    """ Read the table containing the zero points for a given tile and given
    bands. """
    zpfile = os.path.join(context.tables_dir, "ZPfiles_Feb2019",
                          "{}_ZP.cat".format(tile))
    zpdata = Table.read(zpfile, format="ascii")
    zpdict = dict([(t["FILTER"], t["ZP"]) for t in zpdata])
    zps = np.array([zpdict[band] for band in bands])
    return zps
