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
    zpfile = os.path.join(context.tables_dir, "ZPfiles_2020",
                                       "MAIN3.1_ZPs.cat")

    zpdata = Table.read(zpfile, format="ascii")
    zpdic = {a: {"R": b,
                 "F660": c,
                 "I": d}
              for a, b, c, d in zip(zpdata["FIELD"],
                                    zpdata["SPLUS_R"],
                                    zpdata["SPLUS_F660"],
                                    zpdata["SPLUS_I"])} 
    zps = np.array([zpdic[Field][band] for band in bands])
    return zps
