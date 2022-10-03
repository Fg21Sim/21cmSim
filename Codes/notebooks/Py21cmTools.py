# Copyright (c) 2020-2022 Chenxi SHAN <cxshan@hey.com>
# 21cm Tools
__version__ = "0.0.2"
__date__    = "2022-06-12"

import os
import sys
import argparse

from datetime import datetime, timezone
import time

import numpy as np
from scipy import stats

import pickle
from pprint import pprint

from astropy.io import fits
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
import astropy.units as au
from astropy.cosmology import Planck15

import logging
logger = logging.getLogger('21cmFAST')
logger.setLevel(logging.INFO)
import py21cmfast as p21c
# For plotting the cubes, we use the plotting submodule:
from py21cmfast import plotting
# For interacting with the cache
from py21cmfast import cache_tools

#=================== freq & z tools ===================#

def freq21cm( ref_freq, default=True ):
    r"""
    *** Get the rest-frame freq of a 21cm emission in MHz ***
    
    :Params ref_freq: the referenced rest-frame freq of the 21cm line in Hz;
    :Params default: default = False uses the ref_freq, else uses the default value 1420405751.7667 Hz;
    :Output freq21cm: the rest-frame freq of a 21cm emission in MHz;
    """
    if default==False:
        freq21cm = ref_freq / 1e6  # [MHz]
    else:
        freq21cm = 1420405751.7667 / 1e6  # [MHz]
    return freq21cm

def fz_list( begin, step, stop ):
    r"""
       *** Generate a list of freqs or z *** # Basic function
       !!!  !!! # Important note
       +++  +++ # Improvements
       :Params begin: start of the list
       :Params step: step of the list
       :Params stop: end of the list
       :Output values: freqs or z list
   """
    values = []
    begin, step, stop = float(begin), float(step), float(stop)
    v = np.arange(start=begin, stop=stop+step/2, step=step)
    values += list(v)
    return values

def z2freq( redshifts, print_=False ):
    r"""
    *** Convert the redshift to freq ***
    :Params redshifts: input redshift;
    :Params print_: print indicator;
    :Output freqs: output freq [a list of float];
    """
    redshifts = np.asarray(redshifts)
    freqs = freq21cm(1) / (redshifts + 1.0)
    if print_:
        print("# redshift  frequency[MHz]")
        for z, f in zip(redshifts, freqs):
            print("%.4f  %.2f" % (z, f))
    return freqs

def freq2z( freqs, print_=False ):
    r"""
       *** Convert the redshift to freq *** # Basic function
       !!!  !!! # Important note
       +++  +++ # Improvements
       :Params freqs: input freqs;
       :Params print_: print indicator;
       :Output redshifts: output redshifts [a list of float];
   """
    freqs = np.asarray(freqs)
    redshifts = freq21cm(1) / freqs - 1.0
    if print_:
        print("# frequency[MHz]  redshift")
        for f, z in zip(freqs, redshifts):
            print("%.2f  %.4f" % (f, z))
    return redshifts

def wlist( list_, filename ):
    r"""
       *** Write a list to file w/ line break *** # Basic function
       !!!  !!! # Important note
       +++  +++ # Improvements
       :Params list_: python list object;
       :Params filename: path & name to the file;
       :Output None:
   """
    with open(filename, mode="w") as outfile:
        for s in list_[:-1]:
            outfile.write("%s\n" % s)
        outfile.write("%s" % s)
    print(filename, 'is saved!')

def genfz( type_, begin, step, stop, print_=False ):
    r"""
       *** Generate a freq & z dict *** # Basic function
       !!!  !!! # Important note
       +++  +++ # Improvements
       :Params type_: 'z' redshift; 'f' frequency;
       :Params begin: start of the list;
       :Params step: step of the list;
       :Params stop: end of the list;
       :Output fz: freq & z dict;
   """
    fz = {}
    inlist = fz_list( begin, step, stop )
    if type_ == 'z':
        fz['z'] = inlist
        fz['f'] = z2freq( inlist )
    elif type_ == 'f':
        fz['f'] = inlist
        fz['z'] = freq2z( inlist )
    else:
        print("Check your list type!")
    if print_:
        for z, f in zip(fz['z'], fz['f']):
            print("z=%06.3f, freq=%06.2f MHz" % (z, f))
    return fz
    
    
#=================== Cube releated tools ===================#

def swap_axes( lightbox, axes_swapped=True ):
    r"""
       *** Swap axis of the lightbox *** # Basic function
       !!! lightbox = brightness temp  !!! # Important note
       +++  +++ # Improvements
       :Params lightbox: brightness temp
               lightbox shape:
                  * [X, Y, LoS] (lightcone object shape)
                  * [LoS, Y, Y] (lighttravel shape)
       :Params axes_swapped: Only swap axis=0 & axis=2 when True;
       :Output swapped: Only give swapped lightbox when axes_swapped = True;
   """
        
    if axes_swapped == True:
        swapped = np.swapaxes(lightbox.copy(), 0, 2)
    else:
        swapped = lightbox
    return swapped

#=================== Lightcone releated tools ===================#

def loadlightcone( fname='Lightcone200M294', direc='/home/cxshan/DataBank/EoR/Custom21cmSim' ):
    r"""
        *** Load lightcone *** # Basic function
        !!!  !!! # Important note
        +++  +++ # Improvements
        :Params fname: file name
        :Params direc: file dir
        :Output lightcone: 21cmfast Lightcone object
    """
    lightcone = p21c.outputs.LightCone.read( fname=fname,direc=direc )
    return lightcone

def lightconeinfo( lightcone ):
    r"""
    *** Save the lightcone info in a dict ***
    :Params lightcone: the input lightcone;
    :Output info: the info dict contains all the Lightcone attribute;
    """
    info = {}
    info['cell_size'] = lightcone.cell_size
    info['global_xHI'] = lightcone.global_xHI
    info['lightcone_coords'] = lightcone.lightcone_coords
    info['lightcone_dimensions'] = lightcone.lightcone_dimensions
    info['lightcone_distances'] = lightcone.lightcone_distances
    info['n_slices'] = lightcone.n_slices
    info['lightcone_redshifts'] = lightcone.lightcone_redshifts
    info['shape'] = lightcone.shape
    return info

def pklsave( data, file ):
    r"""
    *** Save intermediate files w/ pkl ***
    Params data: data to be stored;
    Params file: file name to be stored;
    """
    #Save skymodel
    pickle_file = open(file, 'wb')
    pickle.dump(data, pickle_file,  protocol=4)
    pickle_file.close()
    print(file, "is saved!")
    
def lightcone_voxels( lightcone, print_=True ):
    r"""
        *** Calculate the voxels *** # Basic function
        !!! Normally the cell & slice size of a voxel is the same; !!! # Important note
        +++ Add a way to check the axes_swapped keywords +++ # Improvements
        :Params lightcone: 21cmFAST lightcone;
        :Params print_: Print indicator;
        :Output voxels: voxel size numpy array in [Side-X,Side-Y,LoS] or [LoS,Side-Y,Side-X]
                        if axes_swapped!
    """
    voxels = np.asarray(lc.lightcone_dimensions) / np.asarray(lc.shape)
    if print_:
        print("The Cell & Slice size [Mpc] of the lightcone voxels:")
        print("- dDc_side=%.3f [cMpc], dDc_slice=%.3f [cMpc]" % (voxels[1], voxels[0]))
    return voxels


#=================== FgSimLightCone tools ===================#

class SimLightCone:
    r"""
       *** Create a FgSimLightCone from a 21cmFAST lightcone *** # Basic function
       !!! In active dev mode !!! # Important note
       +++ Add axis swap +++ # Improvements
       :Params fname:
       :Params direc:
       :Params unit:
       :Params recalculate_:
       :Params default_:
       :Params print_:
    """
    def __init__(self, fname, direc, unit='K', recalculate_=False, default_=True, print_=False):
        self.unit = unit
        self.recalculate_ = recalculate_
        self.default_ = default_
        self.print_ = print_
        self.fname = fname
        self.direc = direc
        # Load the 21cmFAST lightcone
        self.FASTlightcone = p21c.outputs.LightCone.read( fname=fname,direc=direc )
        logger.info("Loaded light-cone cube: %dx%d (cells) * %d (slices)" %
                    (self.Nside, self.Nside, self.Nslice))\

    @property
    def voxels( self ):
        r"""
            *** Calculate the voxels *** # Basic function
            !!! Normally the cell & slice size of a voxel is the same; !!! # Important note
            +++ Add a way to check the axes_swapped keywords +++ # Improvements
            :Params lightcone: 21cmFAST lightcone;
            :Params print_: Print indicator;
            :Output voxels: voxel size numpy array in [Side-X,Side-Y,LoS] or [LoS,Side-Y,Side-X]
                            if axes_swapped!
        """
        voxels = np.asarray(self.FASTlightcone.lightcone_dimensions) / np.asarray(self.FASTlightcone.shape)
        logger.info("The Cell & Slice size [Mpc] of the lightcone voxels: dDc_side=%.3f [cMpc], dDc_slice=%.3f [cMpc]" % (voxels[1], voxels[0]))
        if self.print_:
            print("The Cell & Slice size [Mpc] of the lightcone voxels:")
            print("- dDc_side=%.3f [cMpc], dDc_slice=%.3f [cMpc]" % (voxels[1], voxels[0]))
        return voxels

    @property
    def plank18( self ):
        '''
        # Cosmology is from https://arxiv.org/pdf/1807.06209.pdf
        # Table 2, last column. [TT,TE,EE+lowE+lensing+BAO]
        '''
        Planck18 = Planck15.clone(
            Om0=(0.02242 + 0.11933) / 0.6766**2,
            Ob0=0.02242 / 0.6766**2,
            H0=67.66,
        )
        _defaults_ = {
            "SIGMA_8": 0.8102,
            "hlittle": Planck18.h,
            "OMm": Planck18.Om0,
            "OMb": Planck18.Ob0,
            "POWER_INDEX": 0.9665,}
        return _defaults_

    @property
    def cosmo( self, H0=67.8, Om0=0.308, Ob0=0.0484 ):
        r"""
            *** Calaculate the Cosmology *** # Basic function
            !!! Use the `21cmFAST` lightcone cosmology !!! # Important note
            :Params H0: preset=67 Hubble constant at z = 0.  If a float, must be in [km/sec/Mpc];
            :Params Om0: preset=0.27 Omega matter: density of non-relativistic matter in units of the critical density at z=0.
            :Params Ob0: preset=0.046 Omega baryons: density of baryonic matter in units of the critical density at z=0.
            :Output cosmo: Return an astropy cosmology object.
        """
        if self.default_:
            print("====== Using the Plank18 cosmology ======")
            cosmo = Planck15.clone(H0=self.plank18['hlittle'] * 100, Om0=self.plank18['OMm'], Ob0=self.plank18['OMb'])
            logger.info("Using the Plank18 cosmology: {}".format(cosmo))
        else:
            print("====== Using the Custom cosmology ======")
            cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0)
            logger.info("Using the Custom cosmology: {}".format(cosmo))
        if self.print_:
            print(cosmo)
        return cosmo

    @property
    def Nslice( self ):
        _, _, nslice = self.FASTlightcone.shape
        return nslice

    @property
    def Lslice( self ):
        _, _, lslice = self.FASTlightcone.lightcone_dimensions
        return lslice

    @property
    def Nside( self ):
        _, nside, _ = self.FASTlightcone.shape
        return nside

    @property
    def Lside( self ):
        _, lside, _ = self.FASTlightcone.lightcone_dimensions
        return lside

    @property
    def slices_Dc( self ):
        """
        The comoving distances of each slice in the light-cone cube.
        The slices are evenly distributed along the LoS with equal
        comoving step. [Mpc]
        """
        r"""
            ***  *** # Basic function
            !!!  !!! # Important note
            +++  +++ # Improvements
            :Params recalculate_: recalculate the comoving distance;
            :Output Dc: Comoving distance
        """
        if self.recalculate_:
            print("====== Recalculate the Dc instead of using the original ======")
            dDc_side, _, dDc_slice = self.voxels
            Dc_step = dDc_slice
            Dc_min = self.FASTlightcone.lightcone_distances[0]
            Dc = np.array([Dc_min + Dc_step*i for i in range(self.Nslice)])
            logger.info("Recalculate the Dc instead of using the original.")
            if self.print_:
                print(Dc)
        else:
            Dc = self.FASTlightcone.lightcone_distances
        return Dc

    @property
    def wcs( self ):
        w = WCS(naxis=2)
        w.wcs.ctype = ["pixel", "pixel"]
        w.wcs.cunit = ["Mpc", "Mpc"]  # [cMpc]
        w.wcs.crpix = np.array([1.0, 1.0])
        w.wcs.crval = np.array([0.0, 0.0])
        w.wcs.cdelt = np.array([self.voxels[0], self.voxels[1]]) # [cMpc]
        return w

    @property
    def brightness( self ):
        r"""
        *** Converting brightness temperature default mK to custom unit ***
        :Params self.FASTlightcone.brightness_temp: the input brightness temperature;
        :Params unit: unit to be converted;
        :Output brightness: unit converted brightness temp data;
        """
        default = 'mK'
        unit_in = au.Unit(default)
        unit_out = au.Unit(self.unit)
        brightness = self.FASTlightcone.brightness_temp.copy()
        brightness *= unit_in.to(unit_out)
        if self.unit == 'K':
            unit_converted = True
            self.unit_converted = unit_converted
            logger.info("Brightness is converted to {}".format(self.unit))
        return brightness

    def get_slice( self, z, print_=True ):
        r"""
            *** Get the slice of aimed redshift *** # Basic function
            !!! Assue linear relation for intermediate redshift !!! # Important note
            +++  +++ # Improvements
            :Params z: aimed redshift;
            :Output zslice: slice at given z;
        """
        Dc = self.cosmo.comoving_distance(z).value  # [cMpc]
        if print_:
            print(self.cosmo)
            print(Dc)
        slices_Dc = self.slices_Dc
        if Dc < slices_Dc.min() or Dc > slices_Dc.max():
            raise ValueError("requested redshift out of range: %.2f" % z)

        i2 = (slices_Dc <= Dc).sum()
        i1 = i2 - 1
        Dc1, s1 = slices_Dc[i1], self.brightness[:, :, i1]
        Dc2, s2 = slices_Dc[i2], self.brightness[:, :, i2]
        slope = (s2 - s1) / (Dc2 - Dc1)
        zslice = s1 + slope * (Dc - Dc1)
        if print_:
            print( 'Nearest i:', i2, 'Nearest i-1:', i1 )
            print( 'Dc1:', Dc1, 's1:', s1 )
            print( 'Dc2:', Dc2, 's2:', s2 )
            print( 'slope:', slope, 'zslice:', zslice )
        return zslice

    def write_slice( self, outfile, data, z, freq, clobber=False, print_=True ):
        Dc = self.cosmo.comoving_distance(z).value  # [cMpc]
        # Header info
        header = fits.Header()
        header["BUNIT"] = (self.unit, "Data unit")
        header['FASTlightcone_n'] = (self.fname, "Parent lightcone name")
        header['FASTlightcone_d'] = (self.direc, "Parent lightcone path")
        header['zmin'] = (self.FASTlightcone.lightcone_redshifts[0], "Parent lightcone redshift min")
        header['zmax'] = (self.FASTlightcone.lightcone_redshifts[-1], "Parent lightcone redshift max")
        # Cosmology
        if self.default_:
            header['cosmo'] = ("Default", "Using Plank18 from Plank15 clone w/ updated H0, Om0, & Omb")
        else:
            header['cosmo'] = ("Custom", "Using Custom FlatLambdaCDM w/ H0, Om0, & Omb")
        header['H0'] = (self.cosmo.H0.value, "H0: Hubble constant at z = 0.")
        header['Om0'] = (self.cosmo.Om0, "Omega matter: density of non-relativistic matter in units of the critical density at z=0.")
        header['Ob0'] = (self.cosmo.Ob0, "Omega baryons: density of baryonic matter in units of the critical density at z=0.")
        # Dimensions
        header["Lside"] = (self.Lside,
                           "[cMpc] Simulation cell side length")
        header["Nside"] = (self.Nside,
                           "[N] Number of cell at each cell side")
        # Slice info
        header["REDSHIFT"] = (z, "redshift of this slice")
        header["FREQ"] = (freq, "[MHz] observed HI signal frequency")
        header["Dc"] = (Dc, "[cMpc] comoving distance")
        # Anc
        header["DATE"] = (datetime.now(timezone.utc).astimezone().isoformat(),
                          "File creation date")
        #header.add_history(" ".join(sys.argv))
        # WCS
        header.extend(self.wcs.to_header(), update=True)
        
        if print_:
            print(header)
        hdu = fits.PrimaryHDU(data=data, header=header)
        try:
            hdu.writeto(outfile, overwrite=clobber)
        except TypeError:
            hdu.writeto(outfile, clobber=clobber)
        logger.info("Wrote slice to file: %s" % outfile)