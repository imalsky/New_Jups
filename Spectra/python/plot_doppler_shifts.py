'''

Caleb Harada

plot doppler shiifts as function of phase


'''


import numpy as np
from PyAstronomy import pyasl
from scipy.optimize import curve_fit
import matplotlib.pylab as plt
import glob
import matplotlib.ticker as ticker



### ----- INPUT/OUTPUT ----- ###


savefile = 'fig/doppler_shifts.png'





### ----- GET FILES ----- ###

clear_dopp = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/clear/*_dopp-ON_phase-*.dat'
clear_rest = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/clear/*_dopp-OFF_phase-*.dat'


# TEST CASE (old RT)
#clear_test_dopp = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/clear/*_dopp-ON_phase-*.dat'
#clear_test_rest = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/clear/*_dopp-OFF_phase-*.dat'


compact_thin_dopp = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/compact_thin/*_clouds-ON_dopp-ON_phase-*.dat'
compact_thin_rest = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/compact_thin/*_clouds-ON_dopp-OFF_phase-*.dat'

compact_thick_dopp = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/compact_thick/*_clouds-ON_dopp-ON_phase-*.dat'
compact_thick_rest = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/compact_thick/*_clouds-ON_dopp-OFF_phase-*.dat'

extended_thin_dopp = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/extended_thin/*_dopp-ON_phase-*.dat'
extended_thin_rest = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/extended_thin/*_dopp-OFF_phase-*.dat'

extended_thick_dopp = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/extended_thick/*_dopp-ON_phase-*.dat'
extended_thick_rest = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/extended_thick/*_dopp-OFF_phase-*.dat'

pp_extended_thick_dopp = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/pp_extended_thick/*_dopp-ON_phase-*.dat'
pp_extended_thick_rest = '/Users/calebharada/3D_RT_Code/2stream_new/OUT/pp_extended_thick/*_dopp-OFF_phase-*.dat'





def dopp_signatures(templates, signals, phase_fmt='int'):

	### ----- LOAD DATA ----- ###

	templates = sorted(glob.glob(templates))
	signals = sorted(glob.glob(signals))

	phases = np.zeros(len(templates) + 1)
	RVs = np.zeros(len(templates) + 1)



	# set up CCF figure
	font = {'size' : 14,
        'family' : 'serif'}
	plt.rc('font', **font)

	CCF_fig, ax = plt.subplots(1, 1, figsize=(10, 10))
	ax.set_ylabel('Doppler shift [km s$^{-1}$]')
	ax.set_xlabel('Phase angle [degrees]')
	ax.set_ylim([-10, 10])
	ax.set_xlim([-75, 360])

	ax_cc = ax.twiny()
	ax_cc.set_ylim([-5, 5])
	ax_cc.set_xlim([0, 5.8])
	ax_cc.set_xlabel('Cross correlation + constant')

	# add colorbar
	sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=0, vmax=1))
	sm._A = []
	cbar = plt.colorbar(sm)
	cbar.ax.set_ylabel('Orbital phase')

	# list of colors
	colors = np.linspace(0, 1, len(templates))

	# keep track of color and offset
	color_idx = 0
	offset = 0



	for i in range(len(templates)):


		# save orbital
		if phase_fmt == 'int':
			phase = int(templates[i][-7:-4])
		elif phase_fmt == 'float':
			phase = float(templates[i][-10:-4])
		else:
			print('ERROR\n Unrecognized phase format.')
			exit()

		print(phase)
		phases[i] = phase


		tw, tf = np.loadtxt(templates[i], unpack=True)
		dw, df = np.loadtxt(signals[i], unpack=True)


		'''
		# Plot template and data in real time
		plt.plot(tw, tf, 'b-', label='Rest-frame (template)')
		plt.plot(dw, df, 'r-', label='Shifted')
		plt.legend()
		plt.xlim([2.311e-6, 2.313e-6])
		plt.ylim([0, 5e-9])
		plt.title(f'$\\varphi=${phase}')

		plt.plot()
		plt.show(block=False)
		plt.pause(0.01)
		plt.close()
		'''
		


		# Carry out the cross-correlation.
		# The RV-range is -10 to +10 km/s in steps of 0.1 km/s.
		# The first and last 200 points of the data are skipped.
		rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -10., 10., 0.1, mode='lin', skipedge=200)


		# normalize cc function 
		cc = (cc - min(cc)) / (max(cc) - min(cc))


		# Find the index of maximum cross-correlation function
		maxind = np.argmax(cc)


		# define Gaussian function
		def gaussian(x, a, x0, sigma):
		    return a * np.exp(-(x - x0) * (x - x0) / (2 * sigma * sigma))


		# fit Guassian to peak region of cc function
		rv_range = rv[maxind - 20: maxind + 20]
		cc_range = cc[maxind - 20: maxind + 20]
		popt, pcov = curve_fit(gaussian, rv_range, cc_range, p0=[1, 0, 1])


		# calculate best fit guassian
		rv_fit = np.linspace(rv[maxind - 20], rv[maxind + 20], 1000)
		cc_fit = gaussian(rv_fit, *popt)


		# index of maximum of Gaussian
		max_gauss = np.argmax(cc_fit)


		# save radial velocity
		RVs[i] = rv_fit[max_gauss]


		# plot CCFs
		ax_cc.plot(cc + offset, rv, lw=2, color=plt.cm.jet(colors[color_idx]))
		ax_cc.plot(cc[maxind] + offset, rv[maxind], 'x', color=plt.cm.jet(colors[color_idx]))

		ax_cc.plot(cc_fit + offset, rv_fit, 'k-', alpha=0.5)
		ax_cc.plot(cc_fit[max_gauss] + offset, rv_fit[max_gauss], 'kx', alpha=0.5)

		offset += 0.1
		color_idx += 1


		'''
		# print stuff and plot
		print("Cross-correlation function is maximized at dRV = ", rv_fit[max_gauss], " km/s")
		if rv_fit[max_gauss] > 0.0:
		  print("  A red-shift with respect to the template")
		else:
		  print("  A blue-shift with respect to the template")

		plt.plot(rv, cc, 'b-')
		plt.plot(rv[maxind], cc[maxind], 'bx')

		plt.plot(rv_fit, cc_fit, 'r-')
		plt.plot(rv_fit[max_gauss], cc_fit[max_gauss], 'rx')

		plt.title(f'$\\varphi$ = {phase}')
		plt.xlabel('Doppler shift [km/s]')
		plt.ylabel('Cross correlation')

		plt.show()
		'''
		



	phases[-1] = 360
	RVs[-1] = RVs[0]

	# plot RV line
	ax.set_zorder(2)
	ax.patch.set_visible(False)
	ax.plot(phases, RVs, 'r', lw=2)
	ax.plot(phases - 360, RVs, 'r', lw=2)


	# horizontal line
	ax.axhline(0, color='k', lw=1, ls='--', zorder=0)



	#plt.show()

	#CCF_fig.savefig(savefile, bbox_inches='tight', dpi=300)


	return phases, RVs
	



clear_phases, clear_RVs = dopp_signatures(clear_rest, clear_dopp, phase_fmt='float')

# TEST
#clear_test_phases, clear_test_RVs = dopp_signatures(clear_test_rest, clear_test_dopp, phase_fmt='float')

compact_thin_phases, compact_thin_RVs = dopp_signatures(compact_thin_rest, compact_thin_dopp, phase_fmt='float')
compact_thick_phases, compact_thick_RVs = dopp_signatures(compact_thick_rest, compact_thick_dopp, phase_fmt='float')
extended_thin_phases, extended_thin_RVs = dopp_signatures(extended_thin_rest, extended_thin_dopp, phase_fmt='float')
extended_thick_phases, extended_thick_RVs = dopp_signatures(extended_thick_rest, extended_thick_dopp, phase_fmt='float')
pp_extended_thick_phases, pp_extended_thick_RVs = dopp_signatures(pp_extended_thick_rest, pp_extended_thick_dopp, phase_fmt='float')




### --- MAKE FIGURE ----- ###



font = {'size' : 14,
        'family' : 'serif'}
plt.rc('font', **font)
fig, ax1 = plt.subplots(1, 1, figsize=(12, 6))



# adjust x axis
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(15.0))
ax1.xaxis.set_major_locator(ticker.MultipleLocator(60.0))
ax1.set_xlim([0,360])
ax1.set_xlabel('Phase angle [deg]')

# adjust y axis
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax1.set_ylim([-5.5,5.5])
ax1.set_ylabel('Doppler shift [km s$^{-1}$]')

ax1.tick_params(axis='both',
                which='both',
                direction='in',
                top = True,
                right = True)
ax1.tick_params(axis='both',
                which='major',
                length=5)

# add top axis with orbital phase
ax2 = ax1.twiny()
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax2.tick_params(axis='both',
                which='both',
                direction='in',
                top = True,
                right = True)
ax2.set_xlabel('Orbital phase')


'''
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

'''



ax1.axhline(0, color='black', linestyle=':', lw=1, alpha=0.8)
ax1.plot(clear_phases, clear_RVs, c='gray', lw=3, label='Clear')

# TEST
#ax1.plot(clear_test_phases, clear_test_RVs, '.', c='k', label='Clear (New)')

ax1.plot(compact_thin_phases, compact_thin_RVs, '#a65628', lw=1, label='Compact thin')
ax1.plot(compact_thick_phases, compact_thick_RVs, c='#377eb8', lw=1, label='Compact thick')
ax1.plot(extended_thin_phases, extended_thin_RVs, c='#4daf4a', lw=1, label='Extended thin')
ax1.plot(extended_thick_phases, extended_thick_RVs, c='#f781bf', lw=1, label='Extended thick')
ax1.plot(pp_extended_thick_phases, pp_extended_thick_RVs, c='#f781bf', ls='--', lw=1, label='PP extended thick')





'''
ax1.axvline(180, c='k', lw=1.5, ls='--', zorder=0)
ax1.annotate('eclipse', xy=(180, -3.5), xycoords='data', xytext=(180, -3.5), bbox=dict(boxstyle='square', facecolor='white'), horizontalalignment='center')
'''


ax1.legend(frameon=False, loc=3, ncol=2)

#plt.show()


fig.savefig(savefile, bbox_inches='tight', dpi=300)












