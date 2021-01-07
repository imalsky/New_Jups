
'''

Caleb Harada

Plot all spectra.


'''


### ----- INPUTS/OUTPUTS ----- ###


#cases = ['clear', 'pp_extended_thick', 'extended_thick']
cases = ['extended_thin', 'compact_thick', 'compact_thin']
savefile = 'fig/spectra_all_doppler_2.png'

line_centers = [2.3092694, 2.3096526, 2.3110131, 2.3115399, 2.3119189, 2.3131558]




### ----- IMPORT LIBRARIES ----- ###

import numpy as np
import glob 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from scipy import constants as const
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable




### ----- SET UP FIGURE ----- ###



# colormap
cm_name = 'batlow'
cm_file = np.loadtxt(f'ScientificColourMaps5/{cm_name}/{cm_name}.txt')
my_colors = mcolors.LinearSegmentedColormap.from_list(cm_name, cm_file)

colormap = my_colors


# change font
font = {'size' : 16,
        'family' : 'serif'}
plt.rc('font', **font)



# make figure with axis
fig, axes = plt.subplots(3, 1, figsize=(16,18), sharex=True)

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
  #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.0005))
  #ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))
  ax.xaxis.set_major_locator(ticker.MultipleLocator(0.0010))
  ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0001))


  ### ----- GET FILES ----- ###

  doppler_files = f'/Users/calebharada/3D_RT_Code/2stream_new/OUT/{case}/*_dopp-ON_phase-*.dat'


  # set up color cycle
  color_idx = np.linspace(0, 1, len(glob.glob(doppler_files)))
  const1 = 1


  ### ----- PLOT THE SPECTRA ----- ###

  # plot Doppler-shifted spectra
  for i, file in enumerate(sorted(glob.glob(doppler_files))):

    wl, flux = np.loadtxt(file, unpack=True)  # load data

    # convert units to kW / m^2 / um
    flux = flux * const.c / (wl * wl) * 1e-9

    # convert wavelenths to microns
    wl *= 1e+6

    # normalize
    flux = (flux - min(flux)) / (max(flux) - min(flux))
    flux = flux / np.median(flux)

    ax.plot(wl, flux + const1, color=colormap(color_idx[i]))
    const1 += 0.5

  ax.set_ylim([0.5, 15.9])

  for line in line_centers:
    ax.axvline(line, color='k', ls='--', lw=1, zorder=0, alpha=0.5)

  
  




# add colorbar
sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=0, vmax=360))
sm._A = []
cbar = plt.colorbar(sm, ax=axes.ravel().tolist(), location='right', shrink=0.95, pad=0.02, aspect=40)
cbar.ax.set_ylabel('Orbital phase [deg]')





#axes[0].set_xlim([2.3085,2.3105])
axes[0].set_xlim([2.3085,2.3135])


# set x and y axis lables
axes[2].set_xlabel('Wavelength [$\\mu$m]')
axes[1].set_ylabel('Normalized flux (+ offset)')



'''
# annotate plots
axes[0].annotate('Clear', xy=(2.30905, 9.5),  xycoords='data', 
  xytext=(2.3087, 14.3), horizontalalignment='left', backgroundcolor='w')

axes[1].annotate('Extended Thick (post-proc)', xy=(2.30905, 9.5),  xycoords='data', 
  xytext=(2.3087, 14.3), horizontalalignment='left', backgroundcolor='w')

axes[2].annotate('Extended Thick (active)', xy=(2.30905, 9.5),  xycoords='data', 
  xytext=(2.3087, 14.3), horizontalalignment='left', backgroundcolor='w')
'''


# annotate plots
axes[0].annotate('Extended Thin', xy=(2.30905, 9.5),  xycoords='data', 
  xytext=(2.3087, 14.3), horizontalalignment='left', backgroundcolor='w')

axes[1].annotate('Compact Thick', xy=(2.30905, 9.5),  xycoords='data', 
  xytext=(2.3087, 14.3), horizontalalignment='left', backgroundcolor='w')

axes[2].annotate('Compact Thin', xy=(2.30905, 9.5),  xycoords='data', 
  xytext=(2.3087, 14.3), horizontalalignment='left', backgroundcolor='w')



### ----- SAVE PLOT ----- ###
fig.savefig(savefile, bbox_inches='tight', dpi=300)

#plt.show()


