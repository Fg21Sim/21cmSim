# Copyright (c) 2020-2022 Chenxi SHAN <cxshan@hey.com>
# 21cm Tools
__version__ = "0.0.2"
__date__    = "2022-06-12"

import textwrap
from pprint import pprint

def p_attribute( obj ):
    pprint(dir( obj ))

def p_vars( obj ):
    pprint(vars( obj ))
    
def p_info_keywords():
    print("""
    cell_size:            Cell size [Mpc] of the lightcone voxels.
    global_xHI:           Global neutral fraction function.
    lightcone_coords:     Co-ordinates [Mpc] of each cell along the redshift axis
    lightcone_dimensions: Lightcone size over each dimension -- tuple of (x,y,z) in Mpc.
    lightcone_distances:  Comoving distance to each cell along the redshift axis, from z=0.
    lightcone_redshifts:  Redshift of each cell along the redshift axis.
    n_slices:             Number of redshift slices in the lightcone.
    shape:                Shape of the lightcone as a 3-tuple.
    """)
    
def p_21cm_keywords():
    print(textwrap.dedent("""
    # Notes Start ↓↓↓↓↓↓↓↓↓↓============↓↓↓↓↓↓↓↓↓↓ 
    
    - Side/Cell as a reference to the cell side;
    - Slice/LoS as a reference to the LoS/redshift side;

    # Notes Ended ↑↑↑↑↑↑↑↑↑↑============↑↑↑↑↑↑↑↑↑↑ 

    # Basic Fits Header Keywords ↓↓↓↓↓↓↓↓↓↓============↓↓↓↓↓↓↓↓↓↓

    BUNIT:      Data unit
    DATE:       File creation date
    Author:     Author of the data
    Cubed:      Cubed sturcture [LoS,Side-Y,Side-X] or [Side-X,Side-Y,LoS]

    Lside:      [cMpc] Simulation cell side length
    Nside:      [N] Number of cell at each cell side
    Lslice:     [cMpc] Simulation LoS length
    Nslice:     [N] Number of slice along LoS
    zmin:       [N] Min Redshift
    zmax:       [N] Max Redshift

    dDc_side:   [cMpc] delta comoving distance along cell side
    dDc_slice:  [cMpc] delta comoving distance along LoS
    Dc_zmin:    [cMpc] comoving distance at zmin
    Dc_zmax:    [cMpc] comoving distance at zmax
    Dc_zstep:   [cMpc] comoving distance between slices
    
    ↓↓↓ For 2D slice:
    REDSHIFT:   [N] Redshift of this slice
    FREQ:       [MHz] Observed HI signal frequency
    Dc:         [cMpc] Comoving distance of this slice
    
    # Basic Fits Header Keywords ↑↑↑↑↑↑↑↑↑↑============↑↑↑↑↑↑↑↑↑↑
    
    # More keywords by Chenxi ↓↓↓↓↓↓↓↓↓↓============↓↓↓↓↓↓↓↓↓↓
    
    header['FASTlightcone_n'] = (self.fname, "Parent lightcone name")
    header['FASTlightcone_d'] = (self.direc, "Parent lightcone path")
    header['cosmo'] = ("Default", "Using Plank18 from Plank15 clone w/ updated H0, Om0, & Omb")
    header['cosmo'] = ("Custom", "Using Custom FlatLambdaCDM w/ H0, Om0, & Omb")
    header['H0'] = (self.cosmo.H0.value, "H0: Hubble constant at z = 0.")
    header['Om0'] = (self.cosmo.Om0, "Omega matter: density of non-relativistic matter in units of the critical density at z=0.")
    header['Ob0'] = (self.cosmo.Ob0, "Omega baryons: density of baryonic matter in units of the critical density at z=0.")
    
    # More keywords by Chenxi ↑↑↑↑↑↑↑↑↑↑============↑↑↑↑↑↑↑↑↑↑
    
    # Global Keywords ↓↓↓↓↓↓↓↓↓↓============↓↓↓↓↓↓↓↓↓↓

    - axes_swapped: True or False
    - unit_convert: True or False

    # Global Keywords ↑↑↑↑↑↑↑↑↑↑============↑↑↑↑↑↑↑↑↑↑

    # WCS Start ↓↓↓↓↓↓↓↓↓↓============↓↓↓↓↓↓↓↓↓↓
    
    For 3D cubes:
        - w.wcs.ctype = ["pixel", "pixel", "pixel"]
        - wcs.cunit = ["Mpc", "Mpc", "Mpc"]  # comoving distance
        - w.wcs.crpix = np.array([1.0, 1.0, 1.0])
        - w.wcs.crval:
            - [Side,Side,LoS] = np.array([0.0, 0.0, Dc_min])
            - [LoS,Side,Side] = np.array([Dc_min, 0.0, 0.0])
        - w.wcs.cdelt:
            - [Side,Side,LoS] = np.array([dDc_cell, dDc_cell, dDc_slice])
            - [LoS,Side,Side] = np.array([dDc_slice, dDc_cell, dDc_cell])

            !!! Important note:
            - To have LoS as the 3rd direction in fits, we need to have 
            LoS in 1st direction in np;
            "In Numpy, the order of the axes is reversed so the first 
            dimension in the FITS file appears last, the last dimension 
            appears first and so on."
            > github.com/ChenxiSSS/Coding/issues/101
            > docs.astropy.org/en/stable/visualization/wcsaxes/slicing_datacubes.html
    For 2D images:
        - w = WCS(naxis=2)
        - w.wcs.ctype = ["pixel", "pixel"]
        - w.wcs.cunit = ["Mpc", "Mpc"]  # [cMpc]
        - w.wcs.crpix = np.array([1.0, 1.0])
        - w.wcs.crval = np.array([0.0, 0.0])
        - w.wcs.cdelt = np.array([dDc_cell, dDc_cell]) # [cMpc]

    # WCS Ended ↑↑↑↑↑↑↑↑↑↑============↑↑↑↑↑↑↑↑↑↑
    """))