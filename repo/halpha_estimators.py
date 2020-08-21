# -*- coding: utf-8 -*-
"""

Created on 08/04/18

Author : Carlos Eduardo Barbosa

Methods to extract H-alpha in the SPLUS filters according to methods presented
in Villela-Rojo et al. 2015 (VR+15).

"""
from __future__ import division, print_function

import os

import astropy.units as u
import numpy as np
from scipy.interpolate import interp1d

def halpha_2F_method(f660, fr):
    """Simplest method: using r-band to extimate the continuum (Pascual et al.
    2007). See equation (1) from VR+15 """
    # Numerical integration of equation (2)
    # deltas = deltax(6562.8 * u.AA)
    deltar = 1319.50 * u.AA # Value calculated using deltax != VR+2015
    deltaf660 = 151.06 * u.AA # Value calculated using deltax != VR+2015
    return  deltaf660 * (f660 - fr) / (1 - deltaf660 / deltar)

def halpha_3F_method(f660, fr, fi, test=False):
    """Three bands method from VR+2015 (equation (3)). """
    # Constants used in the calculation
    a = 1.29075002747
    betas = {"I": 1.1294371504165924e-07 / u.AA,
             'F660': 0.006619605316773919 / u.AA,
             "R" : 0.0007578642946167448 / u.AA}
    ############################################################################
    # Constants were calculated according to integrals indicated below
    if test:
        lam = 6562.8 * u.AA
        alphas_test = alphax(lam)
        betas_test = betax(lam)
        a_test = (alphas_test["R"] - alphas_test["I"]) / \
                 (alphas_test["F660"] - alphas_test["I"])
        print("Passing test for a: {0[0]}".format(np.isclose(a_test, a)))
        print("Beta values are similar: ")
        for key, val in betas.iteritems():
            print(key, np.isclose(betas[key].value, betas_test[key].value))
    ############################################################################
    numen = (fr - fi) - a * (f660 - fi)
    denom = - betas["F660"] * a + betas["R"]
    return numen / denom

def deltax(lam):
    """ Numerical integration of equation (2) in VR+2015. This was done for
    testing only, but indicates that the values for SPLUS are slightly
    different than those used by JPLUS..
    """
    filters_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               "../tables/filter_curves-master")
    filenames = sorted([_ for _ in os.listdir(filters_dir) if
                        _.endswith("dat")])
    deltax = {}
    for fname in filenames:
        filtername = fname.split(".")[0].replace("F0", "F").replace("SDSS",
                           "").replace("JAVA", "").upper().strip()
        wave, trans = np.loadtxt(os.path.join(filters_dir, fname)).T
        wave *= u.AA
        curve = interp1d(wave, trans, kind="linear", fill_value=0.,
                         bounds_error=False)
        val = np.trapz(trans * wave, wave) / curve(lam) / lam
        if np.isfinite(val):
            deltax[filtername] = val
    return deltax

def alphax(lam):
    """ Numerical integration of equation (4) in VR+2015 for SPLUS system."""
    filters_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               "../tables/filter_curves-master")
    filenames = sorted([_ for _ in os.listdir(filters_dir) if
                        _.endswith("dat")])
    alphax = {}
    for fname in filenames:
        filtername = fname.split(".")[0].replace("F0", "F").replace("SDSS",
                           "").replace("JAVA", "").upper().strip()
        wave, trans = np.loadtxt(os.path.join(filters_dir, fname)).T
        wave *= u.AA
        curve = interp1d(wave, trans, kind="linear", fill_value=0.,
                         bounds_error=False)
        term1 = np.trapz(trans * wave * wave, wave) / curve(lam) / lam
        term2 = np.trapz(trans * wave, wave) / curve(lam) / lam
        val = term1 / term2
        if np.isfinite(val):
            alphax[filtername] = val
    return alphax

def betax(lam):
    """ Determination of second term of equation (4) using deltax. """
    deltas = deltax(lam)
    betas = {}
    for key, value in deltas.iteritems():
        betas[key] = 1. / value
    return betas

if __name__ == "__main__":
    pass
