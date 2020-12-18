'''


update_cloudreport_v2.py



updated 11/30/2020 by Caleb Harada




Need to apply this code to T_P_3D file before using in RT calculation

** Main part of code is at line 450! **



-Interpolates to uniform altitude grid

-Doubles the longitude grid of Mike's cloudreport files. 

-NEW: increases vertical resolution of grid ('NTAU_new')




NOTE:
	
	Columns of input file must be ordered:
		lat, lon, level,
		altitude(m), pressure(bars), temp(k), 
		EW vel(m/s), NS vel, vert vel,
		tau_lw[1] ,g0_lw[1], pi0_lw[1],
		tau_lw[2] ,g0_lw[2], pi0_lw[2],
		tau_lw[3] ,g0_lw[3], pi0_lw[3],
		tau_lw[4] ,g0_lw[4], pi0_lw[4]



'''



### ----- INPUTS AND OUTPUTS ----- ###



old_file = f'/home/imalsky/Desktop/UPS-PLANETS/UPS-BIG-G-CLEAR.txt'
new_file = f'/home/imalsky/Desktop/UPS-PLANETS/UPS-BIG-G-CLEAR-250.txt'


smoothing = False


NLAT = 48

NLON = 96

NTAU = 50

NPARAMS = 12

NLON_new = 96	# for output file

NTAU_new = 250  # for output file







#########################################
#			END USER INPUTS				#
#########################################





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





### ----- IMPORT LIBRARIES ----- ###


import numpy as np

import matplotlib.pyplot as plt

from scipy import interpolate

from scipy.signal import savgol_filter





### ----- FIND NEW ALTITUDE GRID ----- ###


def altitudes(data):

	# data must be an array of dimensions [NLAT x NLON x NTAU x NPARAMS]


	# intial min and max altitudes

	z_min = 0

	z_max = 0


	# find highest and lowest altitude

	for i in range(NLAT):

		for j in range(NLON):

			if data[i][j][0][3] > z_max:

				z_max = data[i][j][0][3]

			if data[i][j][NTAU - 1][3] > z_min:

				z_min = data[i][j][NTAU - 1][3]


	# define new grid of altitudes

	z_grid = np.flip(np.linspace(z_min, z_max, NTAU_new), axis=0)


	# return grid of new altitudes

	return z_grid







### ----- FUNCTION TO DOUBLE ALL DATA ----- ###


def double_lons(data):


	data2 = np.copy(data)

	data3 = []


	# add 360 to copy of longitudes (does not duplicate 360)

	for row in data2:

	    row[1] += 360.

	    if row[1] == 360.:

	        last_row = np.copy(row)

	        last_row[1] += 360.

	        data3.append(last_row)


	# combine copies of data        

	double_data = np.concatenate((data, data2, data3))


	# sort double_data by latitude

	double_data = double_data[double_data[:,0].argsort()]


	# sort data to match original format (sorry, this is messy and probably not efficient)

	for i in range(NLAT):

	    chunk = double_data[int(i * len(double_data) / NLAT) : int((i + 1) * len(double_data) / NLAT)]

	    chunk = chunk[chunk[:,1].argsort()]

	    for j in range(NLON_new):

	        sub_chunk = chunk[int(j * len(chunk) / NLON_new) : int((j +1 ) * len(chunk) / NLON_new)]

	        sub_chunk = sub_chunk[sub_chunk[:,2].argsort()]

	        chunk[int(j * len(chunk) / (NLON_new)) : int((j + 1) * len(chunk) / (NLON_new))] = sub_chunk

	    double_data[int(i * len(double_data) / NLAT) : int((i + 1) * len(double_data) / NLAT)] = chunk


	# convert bars to pascals

	for i, row in enumerate(double_data):

		double_data[i][4] *= 1e+5


    # return doubled data grid

	return double_data







### ----- LINEAR INTERPOLATION OVER ENTIRE GRID ----- ###


def LInterp_1d(data, data_new, z_new, param_col):

	# data must be an array of dimensions [NLAT x NLON x NTAU x NPARAMS]
	# z_grid is array of length NTAU containing new altitude grid points



	for i in range(NLAT):

		for j in range(NLON):


			# old altitude grid at this lat, lon

			z_old = data[i][j][:,3]



			# parameter values on old altitude grid

			param_old = data[i][j][:,param_col]



			# linear interpolation function created from old altitudes and parameter values

			param_interp = interpolate.interp1d(z_old, param_old, kind="linear", bounds_error=False, fill_value=0)



			# apply interpolation function (values not on new grid set to zero)

			param_new = param_interp(z_new)



			# change parameter values in data array to the new interpolated values

			for k in range(NTAU_new):

				data_new[i][j][k][param_col] = param_new[k]



	return data_new








### ----- LOGARITHMIC INTERPOLATION OVER ENTIRE GRID ----- ###


def LogInterp_1d(data, data_new, z_new, param_col, integrate=False):

	# data must be an array of dimensions [NLAT x NLON x NTAU x NPARAMS]
	# z_grid is array of length NTAU containing new altitude grid points




	if integrate == True:


		# integrate

		for i in range(NLAT):

			for j in range(NLON):

				for k in range(NTAU - 1):

					data[i][j][k+1][param_col] = data[i][j][k+1][param_col] + data[i][j][k][param_col]




		# define new half-step grid to interpolate over

		dz = (np.max(z_new) - np.min(z_new)) / (NTAU_new - 1)

		half_grid = np.zeros(NTAU_new + 1)

		half_grid[0] = z_new[0] + (0.5 * dz)

		for n in range(len(half_grid) - 1):

			half_grid[n+1] = half_grid[n] - dz



		# do log interpolation


		for i in range(NLAT):

			for j in range(NLON):


				# old altitude grid at this lat, lon

				z_old = data[i][j][:,3]



				# parameter values on old altitude grid

				param_old = data[i][j][:,param_col]



				# make sure no values fall below some epsilon (zeros are bad!)

				epsilon = 1e-10

				for n in range(len(param_old)):

					if param_old[n] < epsilon:

						param_old[n] = epsilon



				# log of old parameter

				log_param_old = np.log(param_old)



				# linear interpolation function created from old altitudes and log parameter values

				param_interp = interpolate.interp1d(z_old, log_param_old, kind="linear", fill_value="extrapolate")



				# apply interpolation function

				log_param_new = param_interp(half_grid)



				# exponentiate interpolated valued to get back to original format

				param_new = np.exp(log_param_new)

				

				# discretize back to grid and change parameter values in data array to the new interpolated values

				for k in range(NTAU_new):

					data_new[i][j][k][param_col] = param_new[k+1] - param_new[k]





	else:



		for i in range(NLAT):

			for j in range(NLON):


				# old altitude grid at this lat, lon

				z_old = data[i][j][:,3]



				# parameter values on old altitude grid

				param_old = data[i][j][:,param_col]



				# make sure no values fall below some epsilon (zeros are bad!)

				epsilon = 1e-10

				for n in range(len(param_old)):

					if param_old[n] < epsilon:

						param_old[n] = epsilon



				# log of old parameter

				log_param_old = np.log(param_old)



				# linear interpolation function created from old altitudes and log parameter values

				param_interp = interpolate.interp1d(z_old, log_param_old, kind="linear", bounds_error=False, fill_value=-1e10)



				# apply interpolation function (values not on new grid should go to zero)

				log_param_new = param_interp(z_new)



				# exponentiate interpolated valued to get back to original format

				param_new = np.exp(log_param_new)



				# change parameter values in data array to the new interpolated values

				for k in range(NTAU_new):

					data_new[i][j][k][param_col] = param_new[k]




	return data_new

	




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





############################################################
###														 ###
###	----- MAIN PART OF CODE -- CALL FUNCTIONS HERE ----- ###
###														 ###
############################################################




# load data and reshape data into more convenient dimensions

data = np.loadtxt(old_file, skiprows=5)


data = data.reshape((NLAT, NLON, NTAU, NPARAMS))
data_new = np.zeros((NLAT, NLON, NTAU_new, NPARAMS))


print(data_new.shape)




# get new altitude grid

z_grid = altitudes(data)



# smooth temperatures (to mitigate numerical noise in cloudy models)

if smoothing == True:

	for i in range(NLAT):

		for j in range(NLON):

			temps = data[i][j][:,5]

			temps = savgol_filter(temps, 7, 3)		# (data, window size, polynomial degree)

			data[i][j][:,5] = temps




# interpolate pressures onto new grid (logarithmic)

data_new = LogInterp_1d(data, data_new, z_grid, 4)


# interpolate temperature onto new grid (linear)

data_new = LInterp_1d(data, data_new, z_grid, 5)


# interpolate winds onto new grid (linear)

data_new = LInterp_1d(data, data_new, z_grid, 6)
data_new = LInterp_1d(data, data_new, z_grid, 7)
data_new = LInterp_1d(data, data_new, z_grid, 8)


# interpolate optical depths onto new grid (logarithmic)

data_new = LInterp_1d(data, data_new, z_grid, 9)
data_new = LInterp_1d(data, data_new, z_grid, 10)
data_new = LInterp_1d(data, data_new, z_grid, 11)


# lastly, set all altitude grids equal (to new grid) and add lat, lon, level

for i in range(NLAT):

	for j in range(NLON):

		for k in range(NTAU_new):

			data_new[i][j][k][3] = z_grid[k]

			data_new[i][j][k][2] = k + 1

			data_new[i][j][k][1] = data[i][j][0][1]

			data_new[i][j][k][0] = data[i][j][0][0]





'''

plt.plot(data[5][7][:,3], data[5][7][:,9], '.')


data = LogInterp_1d(data, z_grid, 9, integrate=True)

plt.plot(z_grid, data[5][7][:,9], 'x')


plt.show()

'''



# double all data, then save to new output file

np.savetxt(new_file, data_new.reshape(NLAT * NLON * NTAU_new, NPARAMS),
           fmt='%12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E\t')







