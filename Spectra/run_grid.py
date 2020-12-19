import pandas as pd
import numpy as np
from scipy.interpolate import Rbf
from scipy.interpolate import RectBivariateSpline



def run_all_grid(planet_name, phases,inclinations, system_obliquity, NTAU, NLAT, NLON, file_type):
    print ('Running the regridding')

    def df_to_txt(file, df):
        np.savetxt(file, df.values, fmt='%12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E  %12.4E\t')

    planet_file = '../Planets/' + planet_name + '.txt'

    for phase in phases:
        for inc in inclinations:

            # Sometime w is read in as an object, but it's set to 0 anyway
            df = pd.read_csv(planet_file,
                             delimiter=r"\s+",
                             names=('lat', 'lon', 'level',
                                    'alt','pres','temp',
                                    'temp_std',
                                    'u', 'v', 'w',
                                    'aero_sw_tau_1', 'sw_asym_1', 'sw_pi0_1',
                                    'aero_sw_tau_2', 'sw_asym_2', 'sw_pi0_2',
                                    'aero_sw_tau_3', 'sw_asym_3', 'sw_pi0_3',
                                    'aero_sw_tau_4', 'sw_asym_4', 'sw_pi0_4'),
                                    dtype={'lat':float,
                                    'lon':float,
                                    'level':float,
                                    'alt':float,
                                    'pres':float,
                                    'temp':float,
                                    'temp_std':float,
                                    'u':float,
                                    'v':float,
                                    'aero_sw_tau_1':float,
                                    'sw_asym_1':float,
                                    'sw_pi0_1':float,
                                    'aero_sw_tau_2':float,
                                    'sw_asym_2':float,
                                    'sw_pi0_2':float,
                                    'aero_sw_tau_3':float,
                                    'sw_asym_3':float,
                                    'sw_pi0_3':float,
                                    'aero_sw_tau_4':float,
                                    'sw_asym_4':float,
                                    'sw_pi0_4':float,
                                    'incident_frac': float},
                                    low_memory=False,
                                    index_col='lat')

            df['w'] = 0.0

            df = df.drop(['temp_std'], axis=1)

            df['pres'] = 1e5 * df['pres']

            # Get rid of any place where the longitude is equal to 360
            # Not needed by Eliza Code
            #df = df[(df['lon'] != 360)].reset_index()
            df = df.reset_index()


            def get_incident_flux(df, system_obliquity):
                """Add a new column corresponding to the fraction of intensity"""
                
                # Find the incident fraction depending on latitude and longitude
                df['incident_frac'] = np.sin(df.lat * np.pi / 180.0) * np.sin(system_obliquity) + \
                                      np.cos(df.lat * np.pi / 180.0) * np.cos(system_obliquity) * np.cos(df.lon * np.pi / 180.0)

                # Anywhere not on the circle illuminated set to 0
                #df.incident_frac = df.incident_frac.mask(df.incident_frac < 0.0, 0.0)
                return df

            def phase_rotation(df, phase):
                """ Rotate the planet a certain phase and rollover longitude"""
                # Rotate by the given phase
                df.lon = df.lon + phase

                # Rollover the longitude values if they are too big
                df.lon = df.lon.mask(df.lon >= 360.0, df.lon - 360.0)
                return df


            def wind_rot(df):
                """ Calculate the new u wind speed"""
                u = df.u
                v = df.v
                lon = df.lon
                lat = df.lat

                obs_theta_degree = 180.0 - phase
                obs_theta = obs_theta_degree * (np.pi / 180.0)
                obs_theta = 0

                # Convert to radians
                phi = lat * (np.pi / 180.0)
                theta = lon * (np.pi / 180.0)

                # Get the prime cordinates
                phi_prime = np.arcsin(np.cos(inc) * np.sin(phi) -
                                      np.cos(theta) * np.sin(inc) * np.cos(phi))

                theta_prime = np.arctan2(np.sin(theta) * np.cos(phi),
                                         np.cos(inc) * (np.cos(theta) * np.cos(phi)) + np.sin(inc) * np.sin(phi))

                # Update the lat and lon with the new rotation
                lat_prime = phi_prime * (180.0 / np.pi)
                lon_prime = theta_prime * (180.0 / np.pi)

                # I am not including any w terms in the vector rotation
                v_x = u * np.sin(theta - obs_theta) + v * np.cos(theta - obs_theta) * np.sin(phi)
                v_y = -u * np.cos(theta - obs_theta) + (v * np.sin(theta - obs_theta) * np.sin(phi))
                v_z = v * np.cos(phi)

                # All this vector math is in Deryl's google doc
                v_x_prime = -v_x * np.cos(inc) - (v_z * np.sin(inc))
                v_y_prime = -v_y
                v_z_prime = v_x * np.sin(inc) + (v_z * np.cos(inc))

                u_prime = (-v_x_prime) * np.sin(theta_prime) + v_y_prime * np.cos(theta_prime)
                v_prime = np.sin(phi_prime) * ((-v_x_prime * np.cos(theta_prime)) - v_y_prime * np.sin(theta_prime)) + v_z_prime * np.cos(phi_prime)

                # Make everything mod 360
                if lat_prime < -90.0:
                    lat_prime = 180.0 + lat_prime
                elif lat_prime > 90.0:
                    lat_prime = lat_prime - 180.0

                if lon_prime < 0.0:
                    lon_prime = 360.0 + lon_prime
                elif lon_prime >= 360.0:
                    lon_prime = lon_prime - 360.0
                 
                df.lat = lat_prime
                df.lon = lon_prime
                df.u = u_prime
                df.v = v_prime
                return df

            df = get_incident_flux(df, system_obliquity)
            df = phase_rotation(df, phase)
            df = df.apply(wind_rot, axis=1)

            running_df = pd.DataFrame(columns=['lat', 'lon', 'level',
                                               'alt', 'pres', 'temp', 
                                               'u', 'v', 'w',
                                               'aero_sw_tau_1', 'sw_asym_1', 'sw_pi0_1',
                                               'aero_sw_tau_2', 'sw_asym_2', 'sw_pi0_2',
                                               'aero_sw_tau_3', 'sw_asym_3', 'sw_pi0_3',
                                               'aero_sw_tau_4', 'sw_asym_4', 'sw_pi0_4',
                                               'incident_frac'])


            levels = np.linspace(1, NTAU, NTAU)
            for level in levels:
                print ("On Level", level)
                params = ['pres','temp', 'incident_frac',
                          'u', 'v', 'w',
                          'aero_sw_tau_1', 'sw_asym_1', 'sw_pi0_1',
                          'aero_sw_tau_2', 'sw_asym_2', 'sw_pi0_2',
                          'aero_sw_tau_3', 'sw_asym_3', 'sw_pi0_3',
                          'aero_sw_tau_4', 'sw_asym_4', 'sw_pi0_4']
    
                full_df = df[(df['level'] == level)].reset_index(drop=True)
                x_full_df = np.array(list(full_df.lon))
                y_full_df = np.array(list(full_df.lat))

                sub_df = df[(df['level'] == level)].reset_index(drop=True)
                sub_df = sub_df[(sub_df['temp'] > 100)]

                t1 = np.round(np.linspace(min(full_df.lon), max(full_df.lon), int(NLON)), 3)
                t2 = np.round(np.linspace(min(full_df.lat), max(full_df.lat), NLAT), 3)

                xx, yy = np.meshgrid(t1, t2)
                
                columns = t1.astype(np.str)
                rows = t2.astype(np.str)

                x = np.array(list(sub_df.lon))
                y = np.array(list(sub_df.lat))
                
                # Check whether I can just use the normal x and y
                rbf_full_temps = Rbf(x_full_df, y_full_df, np.array(list(full_df.temp)), epsilon=1, smooth=1, function='linear')
                z_full_temps = rbf_full_temps(xx, yy)
                full_temp = pd.DataFrame(data=z_full_temps, index=rows, columns=columns).stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: "temp"}) 

                param_dataframes = []
                for param in params:
                    rbf = Rbf(x, y, list(sub_df[param]), epsilon=1, smooth=1, function='linear')
                    z = rbf(xx, yy)

                    temp = pd.DataFrame(data=z,index=rows,columns=columns).stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: param}) 
                    param_dataframes.append(temp)
                
                big_df = param_dataframes[0]
                for i in range(len(param_dataframes) - 1):
                    big_df = pd.merge(big_df, param_dataframes[i+1], how='left')
                    
                big_df['level'] = level
                big_df['alt'] = min(sub_df.alt)

                big_df['pres'][full_temp['temp'] < 0] = 0
                big_df['temp'][full_temp['temp'] < 0] = 0

                big_df['pres'][big_df['pres'] < 0] = 0
                big_df['temp'][big_df['temp'] < 0] = 0

                big_df['aero_sw_tau_1'][big_df['aero_sw_tau_1'] < 0] = 0
                big_df['aero_sw_tau_2'][big_df['aero_sw_tau_2'] < 0] = 0
                big_df['aero_sw_tau_3'][big_df['aero_sw_tau_3'] < 0] = 0
                big_df['aero_sw_tau_4'][big_df['aero_sw_tau_4'] < 0] = 0

                big_df['sw_asym_1'][big_df['sw_asym_1'] < 0] = 0
                big_df['sw_asym_2'][big_df['sw_asym_2'] < 0] = 0
                big_df['sw_asym_3'][big_df['sw_asym_3'] < 0] = 0
                big_df['sw_asym_4'][big_df['sw_asym_4'] < 0] = 0
                
                for param in params:
                    big_df[param][full_temp['temp'] < 100] = 0
                
                big_df['lat'] = big_df['lat'].astype(float)
                big_df['lon'] = big_df['lon'].astype(float)
                big_df = big_df.sort_values(by=['lat', 'lon'], axis=0, ascending=[True, True])

                # Make it isothermal
                #big_df['temp'][big_df['temp'] > 100] = 2000

                # Merge the dataframes
                frames = [running_df, big_df]
                running_df = pd.concat(frames, sort=True)

            running_df['lat'] = running_df['lat'].astype(float)
            running_df['lon'] = running_df['lon'].astype(float)
            running_df['level'] = running_df['level'].astype(float)

            # Sort the data tables
            running_df = running_df[['lat', 'lon', 'level',
                                     'alt', 'pres', 'temp',
                                     'u', 'v', 'w',
                                     'aero_sw_tau_1', 'sw_asym_1', 'sw_pi0_1',
                                     'aero_sw_tau_2', 'sw_asym_2', 'sw_pi0_2',
                                     'aero_sw_tau_3', 'sw_asym_3', 'sw_pi0_3',
                                     'aero_sw_tau_4', 'sw_asym_4', 'sw_pi0_4',
                                     'incident_frac']]
            running_df = running_df.sort_values(by=['lat', 'lon', 'level'], axis=0, ascending=[True, True, True])


            # Put this in the spectra folder
            # That is where it needs to be run from
            df_to_txt('../Spectra/DATA/init_' + planet_name + '_phase_{}_inc_{}.txt'.format(phase, inc), running_df)