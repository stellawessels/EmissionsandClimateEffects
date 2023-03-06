import pandas as pd
import numpy as np

data_list = ('IAGOS_timeseries_2019050116041591.txt', 'IAGOS_timeseries_2019043020424591.txt',
             'IAGOS_timeseries_2019043004153591.txt', 'IAGOS_timeseries_2019042914412591.txt',
             'IAGOS_timeseries_2019021216295591.txt', 'IAGOS_timeseries_2019021122212591.txt',
             'IAGOS_timeseries_2019021102051591.txt', 'IAGOS_timeseries_2019021011295591.txt')
file_index = 0

# Calculation saturation pressure
T = 250
a = [-6.0969385e3, 2.12409642e1, -2.711193e-2, 1.673952e-5, 2.433502]
es = np.exp(a[0] * T**(-1) + a[1] + a[2] * T + a[3] * T**2 + a[4] * np.log(T))
print(es)

data = pd.read_csv(data_list[file_index], sep=' ', skiprows=70)
timestep = 10 # seconds
LTO_limit = 914.4 # meters (based on 3000 ft ICAO)

## Data Analysis
data = data.query('baro_alt_AC_val < 7') # remove all lines with measuring errors
data = data.loc[:,('baro_alt_AC','air_press_AC','air_temp_AC','H2O_gas_PC2','ground_speed_AC')] # only select useful data

data = data.assign(distance = data['ground_speed_AC'] * timestep) # append column with distance per timestep
data = data.assign(saturation_pressure = np.exp(a[0] * data['air_temp_AC']**(-1) + a[1] + a[2] * data['air_temp_AC']
                                                + a[3] * data['air_temp_AC']**2 + a[4] * np.log(data['air_temp_AC'])))
data = data.assign(vapour_pressure = data['air_press_AC'] * data['H2O_gas_PC2'] * 1e-6)
data_nonLTO = data.query('baro_alt_AC > 914.4')
data_temp = data.query('air_temp_AC < 235')
data_ISSR = data.query('vapour_pressure > saturation_pressure')
print(data_ISSR)

## Calculations & Results
full_distance = data['distance'].sum() / 1000 # distance in km
distance = data_nonLTO['distance'].sum() / 1000

print('Distance total in km', full_distance)
print('Distance total in km non LTO', distance)

#C Temperature treshold
distance_temp = data_temp['distance'].sum() / 1000

temp_percentage = distance_temp / distance
print(temp_percentage)

#D Schmidt-Appleman Criterion



# print(data.describe(include = 'all'))