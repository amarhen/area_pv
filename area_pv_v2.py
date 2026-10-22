#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 16:01:28 2023

Script to obtain a PV diagram in a FITS cube inside a polygon region in the 
CASA format.

@author: antonio
"""

from astropy.coordinates import SkyCoord
from regions import PixCoord
from regions import Regions

import matplotlib.pyplot as plt
import numpy as np
import os

from astropy import wcs
from scipy import constants
from astropy.io import fits
import matplotlib
import matplotlib.colors as colors

def obtain_pv (file, region, cref, chani = None, chanf = None):
    """
    Function to obtain the PV diagram in a CASA region.

    Parameters
    ----------
    file : string
        Path to the FITS cube to extract the PV diagram.
    region : string
        Path to the CASA region.
    cref : 2-element tuple
        Coordinates of the central point to obtain the PV.
    chani : int, optional
        Initial channel of the cube to obtain the PV. Might save some computing time.
        If not specified, the whole cube is used.
    chanf : int, optional
        Final channel of the cube to obtain the PV. Might save some computing time.
        If not specified, the whole cube is used. Needs 'chani'.
        
    Returns
    -------
    offset : list
        List of distances from the central point. arcsec
    vel : list
        List of velocity values. km/s
    intensity: numpy array
        Intensity of the pixels used for the PV. Same units as FITS file.
    pa : list
        Position angle of the points in the PV. Not used currently.

    """
    
    # Get header
    hdulist = fits.open(file)
    header = hdulist[0].header
    w = wcs.WCS(hdulist[0].header, naxis=2)
    
    naxis1 = header['NAXIS1']  # Number of values in each axis
    naxis2 = header['NAXIS2']  # naxis1 and naxis2 should be equal
    naxis3 = header['NAXIS3']
    
    crval1 = header['CRVAL1']
    crval2 = header['CRVAL2']    
    crval3 = header['CRVAL3'] # starting frequency of the cube in Hz
    
    cdelt1 = header['CDELT1']  # Difference between two pixel in deg
    cdelt2 = header['CDELT2']  
    cdelt3 = header['CDELT3']  # Difference between channels of the cube in Hz
    restfreq = header['RESTFRQ'] # Rest frequency (Hz)
    
    # Build coordinate axes
    ra_px = np.arange(start = 0,    # RA axis in pixel
                 stop = naxis1, 
                 step = 1)
    
    dec_px = np.arange(start = 0,   # DEC axis in pixel
                 stop = naxis2,
                 step = 1)
    
    if chani == None:
        freq_hz = np.arange(start = crval3, # FREQ axis in Hz
                        stop = crval3 + naxis3*cdelt3, 
                        step = cdelt3)
    else:
        freq_hz = np.arange(start = crval3 + chani*cdelt3, # FREQ axis in Hz
                        stop = crval3 + chanf*cdelt3, 
                        step = cdelt3)

    vel = [(((restfreq - i))/restfreq)*constants.c/1000 for i in freq_hz] # vel axis in km/s
    
    # Build (RA,DEC) matrix
    rai_px, deci_px = np.meshgrid(ra_px, dec_px)
    
    # Read region from file
    
    # First, make a temporary copy to solve any possible format issue. Leave text only until ICRS
    with open (region, 'r') as i, open ('tempregionfile.txt', 'w') as o:
        for line in i:
            if 'ICRS' in line:
                o.write(line.partition('ICRS')[0]+line.partition('ICRS')[1])
                break
            else:
                o.write(line)
    
    sky_regions = Regions.read('tempregionfile.txt', format='crtf')
    sky_region = sky_regions[0]
    os.remove('tempregionfile.txt') # Erase the temporary file
    
    # Convert to pixel
    pixel_region = sky_region.to_pixel(w)
    
    # Instance of SkyCoord with the reference coordinates
    cref_sky = SkyCoord(cref[0],cref[1], unit='deg', frame='icrs')
    
    # Store offset, frequency and flux of the points inside the region
    offset = []
    pa = []
    intensity = []
    
    for i in range(len(ra_px)): # Loop over the RA axis
        for j in range(len(dec_px)): # Loop over the DEC axis
            if PixCoord(i,j) in pixel_region: # Check if point is inside the region
                
                intensity.append([])
                
                c = PixCoord(i,j).to_sky(w) # Convert pixel coordinate of the point to sky
                sep = c.separation(cref_sky) # Compute distance in deg to reference sky coordinate
                offset.append(sep.value*3600)
                pa.append(c.position_angle(cref_sky))
                
                if chani == None:
                    for k in range(len(freq_hz)): # Loop over the frequency axis
                        # Store separation (arcsec), frequency (Hz) and flux density (Jy(/beam) of each channel for the point
                        intensity[-1].append(hdulist[0].data[k-1][j][i]) # HDUlist has this weird way of indexing RA,DEC pixel coordinates.
                else:
                    for k in range(len(freq_hz)): # Loop over the frequency axis
                        # Store separation (arcsec), frequency (Hz) and flux density (Jy(/beam) of each channel for the point
                        intensity[-1].append(hdulist[0].data[k+chani-1][j][i]) # HDUlist has this weird way of indexing RA,DEC pixel coordinates.
                        
    return offset, vel, np.array(intensity), pa


def plot_pv(offset, vel, intensity, vsys = 0.0, velmin = 0.0, velmax = 1.0, logscale=True, savefig=False, savefilepath=''):
    """
    Function to plot the position-velocity diagram obtained in 'obtain_pv'.

    Parameters
    ----------
    offset : list
        List of distances from the central point. arcsec
    vel : list
        List of velocity values. km/s
    intensity: numpy array
        Intensity of the pixels used for the PV. Same units as FITS file.
    vsys : float, optional
        Systemic velocity of the source. The default is 0.0.
    velmin : float, optional
        Minimum velocity in the plot. The default is 0.0.
    velmax : float, optional
        Maximum velocity in the plot. The default is 1.0.
    logscale : bool, optional
        If True, color logscale is used. Else, linear. The default is True.
    savefig : bool, optional
        If true, the PV plot is stored in memory. The default is False.
    savefilepath : string, optional
        If savefig==True, this is the path to save the figure. The default is ''.

    Returns
    -------
    None.

    """
    
    #TODO average neighbor pixels to reduce computing time (?)
    
    # Build 2D axes of the variables for the scatter plot
    veli, offseti = np.meshgrid(vel, offset)
    
    xmin = min(offset)
    xmax = max(offset)

    plt.figure()
    plt.xlim(xmin, xmax)
    plt.ylim(velmin, velmax)
    
    xaxis_font = {'size':16, 'weight':'normal'}
    yaxis_font = {'size':16, 'weight':'normal'}

    plt.xlabel(r'Offset from central source (arcsec)', **xaxis_font)
    plt.ylabel(r'Velocity (km s$^{-1}$)', **yaxis_font)
    matplotlib.rc('xtick', labelsize=16)
    matplotlib.rc('ytick', labelsize=16)
    plt.tick_params(labelsize=16) 
    
    plt.axhline(vsys, linestyle = '--', color = 'black', linewidth = 2.0) # Plot a dashed line at vsys
    
    arr = np.linspace(0, 50, 100).reshape((10, 10))
    custom_cmap = truncate_colormap(plt.get_cmap('turbo'), 0.15, 1.0)
    
    if logscale:
        plt.scatter(offseti, veli, c=intensity, cmap = custom_cmap, norm=matplotlib.colors.LogNorm(vmin=0.001, vmax=np.nanmax(intensity)), marker = 's', s=75)
    else:
        plt.scatter(offseti, veli, c=intensity, cmap = custom_cmap, marker = 's', s=55)
    
    if savefig:
        plt.savefig(savefilepath, dpi=500, bbox_inches='tight')
        plt.show()
        plt.close()
    else:
        plt.show()
        
        
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
        
#%%
# Path to the FITS file:
# file = '/Users/antonio/Documents/BHB07_FAUST/data/12m+JOINT-TAKESHI/BHB07-11-Setup2-cube-products_validated_newrestfreqs/images/BHB07-11_260.255GHz_H13COp_joint_0.5_clean.image.pbcor.common.fits'
# file = '/Users/antonio/Documents/BHB07_FAUST/data/12m+JOINT-TAKESHI/BHB07-11-Setup2-cube-products_validated_newrestfreqs/images/BHB07-11_262.004GHz_CCH_joint_0.5_clean.image.pbcor.common.fits'
# file = '/Users/antonio/Documents/BHB07_FAUST/data/CO(3-2)/member.uid___A001_X1465_X80d.BHB07-11_sci.spw33.cube.I.pbcor.dropdeg.fits'
file = '/Users/antonio/Documents/BHB07_FAUST/data/12m+JOINT-TAKESHI/BHB07-11-Setup1-cube-products_for_QA/images/BHB07-11_217.822GHz_c-C3H2_joint_0.5_clean.image.pbcor.common.fits'
# Path to the region file:
# region = '/Users/antonio/Documents/BHB07_FAUST/data/12m+JOINT-TAKESHI/casa_regions/polygon_red_fromh13copm0_version5.crtf'
region = '/Users/antonio/Documents/BHB07_FAUST/data/12m+JOINT-TAKESHI/casa_regions/polygon_blue_fromh13copm0_version2.crtf'


cref = (257.8462667,-27.40914167) # reference center for the PV diagram
vsys = 3.6 # systemic velocity of the source, km/s

# Compute the PV diagram
offset, vel, intensity, pa = obtain_pv(file, region, cref)

# Plot the PV diagram
plot_pv(offset, vel, intensity, vsys, velmin = 1.1, velmax = 4.9, logscale=True, savefig=True, savefilepath='/Users/antonio/Documents/papers_own/faust_bhb11-07_outflow/figures/pvs_outflow/PYTHON_FIGURES/c3h2_blueshifted_jet_newrestfreq_customcolorscale.png')