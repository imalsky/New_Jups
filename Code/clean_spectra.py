
'''

Caleb Harada

Trim edges and interpolate to get rid of Nans. 


'''


### ----- INPUTS/OUTPUTS ----- ###


cases = ['clear_test']
n_trim = 30   # number of points to trim from either end of spectrum



### ----- IMPORT LIBRARIES ----- ###

import numpy as np
import glob 
import os
from scipy import constants as const



def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.
    https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]



# adjust tick marks on both axes
for i, case in enumerate(cases):


  ### ----- GET FILES ----- ###
  files = f'/Users/calebharada/3D_RT_Code/stream_new/OUT/{case}/*_phase-*.dat'


  ### ----- CLEAN UP FILES ----- ###
  for i, file in enumerate(sorted(glob.glob(files))):

    # load data
    wl, flux = np.loadtxt(file, unpack=True)

    # trim edges
    flux_good = flux[n_trim:-n_trim]

    # interpolate over Nans
    nans, x = nan_helper(flux_good)
    flux_good[nans]= np.interp(x(nans), x(~nans), flux_good[~nans])

    # make changes to flux array
    flux[:n_trim] = np.median(flux_good)
    flux[-n_trim:] = np.median(flux_good)
    flux[n_trim:-n_trim] = flux_good

    # convert units to kW / m^2 / um
    flux_lamb = flux * const.c / (wl * wl) * 1e-9

    # convert wavelenths to microns
    wl *= 1e+6

    # save new files
    file_name = os.path.basename(file)
    new_dir = f'/Users/calebharada/3D_RT_Code/rt_emission_aerosols_2stream/OUT/CLEANED/{case}'
    if not os.path.exists(new_dir):
      os.makedirs(new_dir)

    new_file = f'{new_dir}/{file_name}'
    
    np.savetxt(new_file, np.transpose([wl, flux_lamb]), header='wl [microns] \t flux [kW / m^2 / micron]')


