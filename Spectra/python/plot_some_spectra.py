'''

Caleb Harada

Main spectra comparison figure for paper


'''

### ----- INPUTS/OUTPUTS ----- ###


cases = ['clear', 'pp_extended_thick', 'extended_thick']
savefile = 'fig/spectra_compare.png'


labels = ['$\\varphi=45^\\circ$', '$\\varphi=135^\\circ$', '$\\varphi=225^\\circ$', '$\\varphi=315^\\circ$']
line_centers = [2.3092694, 2.3096526]




### ----- IMPORT LIBRARIES ----- ###

import numpy as np
import glob 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import constants as const
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable




### ----- SET UP FIGURE ----- ###



#colors = ['firebrick', 'goldenrod', 'seagreen', 'steelblue']
colors = ['#4daf4a', '#f781bf', '#a65628', '#984ea3']


# change font
font = {'size' : 14,
        'family' : 'serif'}
plt.rc('font', **font)



# make figure with axis
fig, axes = plt.subplots(3, 1, figsize=(14,10), sharex=True)

plt.subplots_adjust(wspace=0, hspace=0)




# adjust tick marks on both axes
for i, ax in enumerate(axes):

  case = cases[i]

  ax.tick_params(axis='both',
                which='both',
                direction='in',
                top = True,
                right = True)
  ax.tick_params(axis='both',
                which='major',
                length=5)

  # set spacing of tick marks
  ax.xaxis.set_major_locator(ticker.MultipleLocator(0.0005))
  ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))


  ### ----- GET FILES ----- ###

  doppler_files = f'/Users/calebharada/3D_RT_Code/2stream_new/OUT/{case}/for_paper/*_dopp-ON_phase-*.dat'
  rest_files = f'/Users/calebharada/3D_RT_Code/2stream_new/OUT/{case}/for_paper/*_dopp-OFF_phase-*.dat'


  # set up color cycle
  color_idx = [0, 1, 2, 3]


  ### ----- PLOT THE SPECTRA ----- ###

  # plot Doppler-shifted spectra
  for i, file in enumerate(sorted(glob.glob(doppler_files))):

    wl, flux = np.loadtxt(file, unpack=True)  # load data

    # convert units to kW / m^2 / um
    flux = flux * const.c / (wl * wl) * 1e-9

    # convert wavelenths to microns
    wl *= 1e+6

    ax.plot(wl, flux, lw=2, color='w')
    ax.plot(wl, flux, lw=1, color=colors[color_idx[i]], label=labels[i])




  # reset color cycle
  plt.gca().set_prop_cycle(None)

  # plot rest-frame spectra
  for i, file in enumerate(sorted(glob.glob(rest_files))):

    wl, flux = np.loadtxt(file, unpack=True)  # load data

    # convert units to kW / m^2 / um
    flux = flux * const.c / (wl * wl) * 1e-9

    # convert wavelenths to microns
    wl *= 1e+6

    ax.plot(wl, flux, '-', lw=1, color=colors[color_idx[i]], alpha=0.3, zorder=1)



  for line in line_centers:
    ax.axvline(line, color='k', ls='--', lw=1, zorder=0, alpha=0.5)


  # y axis limits
  ax.set_ylim([10, 420])




axes[2].legend(loc=1, ncol=2, frameon=False)



# set x limits
axes[2].set_xlim([2.3085,2.3105])


# set x and y axis lables
axes[2].set_xlabel('Wavelength [$\\mu$m]')
axes[1].set_ylabel('Planet flux [kW m$^{-2}$ $\\mu$m$^{-1}$]')


# annotate plots
axes[0].annotate('Clear', xy=(2.30854, 350),  xycoords='data', 
  xytext=(2.30854, 380), horizontalalignment='left')

axes[1].annotate('Post-processed Clouds', xy=(2.30854, 350),  xycoords='data', 
  xytext=(2.30854, 380), horizontalalignment='left')

axes[2].annotate('Active Clouds', xy=(2.30854, 350),  xycoords='data', 
  xytext=(2.30854, 380), horizontalalignment='left')





### ----- SAVE PLOT ----- ###
fig.savefig(savefile, bbox_inches='tight', dpi=300)




