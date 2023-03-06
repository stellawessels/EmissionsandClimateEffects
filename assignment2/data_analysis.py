import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_list = ('IAGOS_timeseries_2019050116041591.txt', 'IAGOS_timeseries_2019043020424591.txt',
             'IAGOS_timeseries_2019043004153591.txt', 'IAGOS_timeseries_2019042914412591.txt',
             'IAGOS_timeseries_2019021216295591.txt', 'IAGOS_timeseries_2019021122212591.txt',
             'IAGOS_timeseries_2019021102051591.txt', 'IAGOS_timeseries_2019021011295591.txt')
file_index = 0

a = [-6.0245282e3, 2.932707e1, 1.0613868e-2, -1.3198825e-5, -4.9382577e-1]

data = pd.read_csv(data_list[file_index], sep=' ', skiprows=70)
timestep = 10 # seconds

## Data Analysis
# data = data.query('baro_alt_AC_val < 7') # remove all lines with measuring errors
data = data.loc[:,('air_press_AC','air_temp_AC','H2O_gas_PC2','ground_speed_AC')] # only select useful data

data = data.assign(distance = data['ground_speed_AC'] * timestep) # append column with distance per timestep
data = data.assign(temperature = data['air_temp_AC'] - 273)

data_nonLTO = data.query('air_press_AC < 75000')
data_temp = data.query('air_temp_AC < 235')

# data = data.query('H2O_gas_PC2 < 0') # remove all lines with measuring errors
data = data.assign(saturation_pressure = np.exp(a[0] * data['air_temp_AC']**(-1) + a[1] + a[2] * data['air_temp_AC']
                                                + a[3] * data['air_temp_AC']**2 + a[4] * np.log(data['air_temp_AC'])))
data = data.assign(vapour_pressure = data['air_press_AC'] * data['H2O_gas_PC2'] / 1e6)

data_ISSR = data.query('vapour_pressure > saturation_pressure')

## Calculations & Results
distance = data['distance'].sum() / 1000
nonLTO_distance = data_nonLTO['distance'].sum() / 1000

distance_temp = data_temp['distance'].sum() / 1000
distance_ISSR = data_ISSR['distance'].sum() / 1000

temp_percentage = distance_temp / distance * 100
ISSR_percentage = distance_ISSR / distance * 100

print(temp_percentage, "percentage of distance under temperature constraint")
print(ISSR_percentage, "percentage of distance under ISSR constraint")

print(data_ISSR.describe(include = 'all'))