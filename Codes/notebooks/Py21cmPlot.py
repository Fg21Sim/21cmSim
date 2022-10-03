# Copyright (c) 2020-2022 Chenxi SHAN <cxshan@hey.com>
# 21cm Tools
__version__ = "0.0.2"
__date__    = "2022-06-12"

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import colors

from py21cmfast import plotting


def EoR_color():
    r"""
       *** Add EoR cmap from 21cmFast  *** # Basic function
       !!! No need to use it if you have 21cmFast !!! # Important note
   """
    eor_colour = colors.LinearSegmentedColormap.from_list(
        "EoR",
        [
            (0, "white"),
            (0.21, "yellow"),
            (0.42, "orange"),
            (0.63, "red"),
            (0.86, "black"),
            (0.9, "blue"),
            (1, "cyan"),
        ],
    )
    plt.register_cmap(cmap=eor_colour)


def plot_Tb(array, unit="mK", label="Brightness Temp", maxrange = 30, minrange = -150):
    r"""
       *** Plot Lightcone slice like p21c.plotting.lightcone_sliceplot *** # Basic function
       !!! max & min range are taken from p21c.plotting.lightcone_sliceplot !!! # Important note
       +++ Add z ticks +++ # Improvements
       :Params array: input slices not lightcone!
       :Params unit: using 'mK' as default, use 'K' for fg21sim; 
       :Params label: label;
       :Params maxrange: max of the colormap;
       :Params minrange: min of the colormap;
   """
    fig = plt.figure(dpi=300)
    sub_fig = fig.add_subplot(111)
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=10)
    c_dens = sub_fig.imshow(array,cmap="EoR")
    
    if unit=='mK':
        minrange = minrange/1
        maxrange = maxrange/1
        c_dens.set_clim(vmin=minrange,vmax=maxrange)
        c_bar = fig.colorbar(c_dens, orientation='vertical')
        c_bar.set_label(r'${\rm \delta T_b [\mathrm{mK}]}$', fontsize=10, rotation=-90, labelpad=10)
    elif unit=='K':
        minrange = minrange/1000
        maxrange = maxrange/1000
        c_dens.set_clim(vmin=minrange,vmax=maxrange)
        c_bar = fig.colorbar(c_dens, orientation='vertical')
        c_bar.set_label(r'${\rm \delta T_b [\mathrm{K}]}$', fontsize=10, rotation=-90, labelpad=10)

    tick_array = np.linspace(minrange, maxrange, 8)
    c_bar.set_ticks(tick_array)


def sview( array ):
    r"""
        *** Imshow an array w/ colorbar  *** # Basic function
        !!!  !!! # Important note
        +++  +++ # Improvements
        :Params array: array to plot;
        :Params :
        :Params :
        :Params :
    """
    plt.imshow(array)
    plt.colorbar()