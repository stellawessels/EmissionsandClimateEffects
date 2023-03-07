import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_list = ('IAGOS_timeseries_2019050116041591.txt', 'IAGOS_timeseries_2019043020424591.txt',
             'IAGOS_timeseries_2019043004153591.txt', 'IAGOS_timeseries_2019042914412591.txt',
             'IAGOS_timeseries_2019021216295591.txt', 'IAGOS_timeseries_2019021122212591.txt',
             'IAGOS_timeseries_2019021102051591.txt', 'IAGOS_timeseries_2019021011295591.txt')
file_index = 3

# Calculation saturation pressure
a = [-6.0245282e3, 2.932707e1, 1.0613868e-2, -1.3198825e-5, -4.9382577e-1]

# Parameters
cp = 1.0035e3 #J/kg/K
M_air = 28.97 # kg/kmol
M_H20 = 18 # kg/kmol
n_1 = 0.3
n_2 = 0.5
Q_1 = 43.2e6 # J/kg
Q_2 = 101.0e6 # J/kg
EI_1 = 1.25 # kg/kg
EI_2 = 9 # kg/kg

data = pd.read_csv(data_list[file_index], sep=' ', skiprows=70)
timestep = 10 # seconds

## Data Analysis
data = data.query('baro_alt_AC_val < 7') # remove all lines with measuring errors
data = data.loc[:,('baro_alt_AC','air_press_AC','air_temp_AC','H2O_gas_PC2','ground_speed_AC')] # only select useful data

data = data.assign(distance = data['ground_speed_AC'] * timestep) # append column with distance per timestep
data = data.assign(saturation_pressure = np.exp(a[0] * data['air_temp_AC']**(-1) + a[1] + a[2] * data['air_temp_AC']
                                                + a[3] * data['air_temp_AC']**2 + a[4] * np.log(data['air_temp_AC'])))
data = data.assign(vapour_pressure = data['air_press_AC'] * data['H2O_gas_PC2'] / 10**6 )
data_nonLTO = data.query('air_press_AC <= 750e2')
data_temp = data.query('air_temp_AC < 235')
data_ISSR = data.query('vapour_pressure > saturation_pressure')

# Calculation liquid 'over water' saturation pressure
a_liq = [-6.0969385e3, 2.12409642e1, -2.711193e-2, 1.673952e-5, 2.433502]
data = data.assign(e_liq = np.exp(a_liq[0] * data['air_temp_AC']**(-1) + a_liq[1] + a_liq[2] * data['air_temp_AC'] +
                                  a_liq[3] * data['air_temp_AC']**2 + a_liq[4] * np.log(data['air_temp_AC'])))

print(data)
## Calculations & Results

#A Flown distance
full_distance = data['distance'].sum() / 1000 # distance in km
distance = data_nonLTO['distance'].sum() / 1000

print('A) Distance total in km', full_distance)
print('Distance total in km non LTO', distance)

#B ISSR

distance_ISSR = data_ISSR['distance'].sum() / 1000
ISSR_percentage = distance_ISSR / distance * 100
print('B) ISSR conditions = ', ISSR_percentage, '%')

#C Temperature treshold
distance_temp = data_temp['distance'].sum() / 1000
temp_percentage = distance_temp / distance * 100
print('C) Temperature Threshold = ', temp_percentage, '%')

#D Schmidt-Appleman Criterion aircraft 1

data = data.assign(G_1= data['air_press_AC'] * cp * (M_air/M_H20) * (EI_1/((1-n_1)*Q_1)))
data = data.assign(T_LM_1 = -46.46 + 9.43* np.log(data['G_1']-0.053) + 0.720 * (np.log(data['G_1'] - 0.053))**2 + 273)
data = data.assign(e_liq1 = np.exp(a_liq[0] * data['T_LM_1']**(-1) + a_liq[1] + a_liq[2] * data['T_LM_1'] +
                                  a_liq[3] * data['T_LM_1']**2 + a_liq[4] * np.log(data['T_LM_1'])))
data = data.assign(LeftSAC_1= data['vapour_pressure'] + data['G_1'] * (data['T_LM_1']-data['air_temp_AC']))
data_SAC_1 = data.query('LeftSAC_1 > e_liq1')

data = data.assign(G_2= data['air_press_AC'] * cp * (M_air/M_H20) * (EI_2/((1-n_2)*Q_2)))
data = data.assign(T_LM_2 = -46.46 + 9.43* np.log(data['G_2']-0.053) + 0.720 * (np.log(data['G_2'] - 0.053))**2 + 273)
data = data.assign(e_liq2 = np.exp(a_liq[0] * data['T_LM_2']**(-1) + a_liq[1] + a_liq[2] * data['T_LM_2'] +
                                  a_liq[3] * data['T_LM_2']**2 + a_liq[4] * np.log(data['T_LM_2'])))
data = data.assign(LeftSAC_2= data['vapour_pressure'] + data['G_2'] * (data['T_LM_2']-data['air_temp_AC']))
data_SAC_2 = data.query('LeftSAC_2 > e_liq2')

"""
data= data.query('air_temp_AC < 235')
plt.plot(data["air_temp_AC"], data["saturation_pressure"])
plt.plot(data["air_temp_AC"], data["vapour_pressure"])
plt.show()
"""
distance_SAC_1 = data_SAC_1['distance'].sum() / 1000
distance_SAC_1_percentage = distance_SAC_1 / distance * 100
print('D) Distance SAC for 1 = ', distance_SAC_1_percentage, '%')

distance_SAC_2 = data_SAC_2['distance'].sum() / 1000
distance_SAC_2_percentage = distance_SAC_2 / distance * 100
print('Distance SAC for 2 = ', distance_SAC_2_percentage, '%')

#E Contrail formation
merged_df_notfinished_1 = pd.merge(data_SAC_1, data_ISSR)
merged_df_1 = pd.merge(merged_df_notfinished_1, data_temp)

merged_df_notfinished_2 = pd.merge(data_SAC_2, data_ISSR)
merged_df_2 = pd.merge(merged_df_notfinished_2, data_temp)

distance_merged_df_1 = merged_df_1['distance'].sum()/1000
distance_merged_df_1_percentage = distance_merged_df_1 / distance * 100

distance_merged_df_2 = merged_df_2['distance'].sum()/1000
distance_merged_df_2_percentage = distance_merged_df_2 / distance * 100
print('E) Persistent contrails for 1 = ', distance_merged_df_1_percentage, '%')
print('Persistent contrails for 2 = ', distance_merged_df_2_percentage, '%')

#F

Bcondition = []

for i in range(0, data.shape[0]):
    if data.iloc[i,7] > data.iloc[i,6]:
        Bcondition.append(1)
    else:
        Bcondition.append(0)

Ccondition =[]

for i in range(0, data.shape[0]):
    if data.iloc[i, 2] < 235:
        Ccondition.append(1)
    else:
        Ccondition.append(0)

Dcondition1 =[]

for i in range(0, data.shape[0]):
    if data.iloc[i, 12] > data.iloc[i, 11]:
        Dcondition1.append(1)
    else:
        Dcondition1.append(0)

Dcondition2 =[]

for i in range(0, data.shape[0]):
    if data.iloc[i, 16] > data.iloc[i, 15]:
        Dcondition2.append(1)
    else:
        Dcondition2.append(0)

# For aircraft 1
ifconditions_1 = []
for i in range(0, data.shape[0]):
    if Bcondition[i] == Ccondition[i] == Dcondition1[i]:
        ifconditions_1.append('Conditions fulfilled')
    else:
        ifconditions_1.append('X')

fault_1 = ['X']

for i in range (1, data.shape[0]-1):
    if ifconditions_1[i] == 'Conditions fulfilled':
        if ifconditions_1[i+1] == 'X':
            fault_1.append('Change')
        else:
            fault_1.append('X')
    else:
        fault_1.append('X')
fault_1.append('X')

wiewasdemol_1 = []

for i in range(0, data.shape[0]):
    if fault_1[i] == 'Change':
        if Bcondition[i+1] == 0:
            wiewasdemol_1.append('B!')
        if Ccondition[i+1] == 0:
            wiewasdemol_1.append('C!')
        if Dcondition1[i+1] == 0:
            wiewasdemol_1.append('D1!')

num_B1s = wiewasdemol_1.count('B!')
num_C1s = wiewasdemol_1.count('C!')
num_D1s = wiewasdemol_1.count('D1!')

sum_1 = num_B1s + num_C1s + num_D1s

percentage_B_1 = num_B1s / sum_1 * 100
percentage_C_1 = num_C1s / sum_1 * 100
percentage_D_1 = num_D1s / sum_1 * 100

# For aircraft 2
ifconditions_2 = []
for i in range(0, data.shape[0]):
    if Bcondition[i] == Ccondition[i] == Dcondition2[i]:
        ifconditions_2.append('Conditions fulfilled')
    else:
        ifconditions_2.append('X')

fault_2 = ['X']

for i in range (1, data.shape[0]-1):
    if ifconditions_2[i] == 'Conditions fulfilled':
        if ifconditions_2[i+1] == 'X':
            fault_2.append('Change')
        else:
            fault_2.append('X')
    else:
        fault_2.append('X')
fault_2.append('X')

wiewasdemol_2 = []

for i in range(0, data.shape[0]):
    if fault_2[i] == 'Change':
        if Bcondition[i+1] == 0:
            wiewasdemol_2.append('B!')
        if Ccondition[i+1] == 0:
            wiewasdemol_2.append('C!')
        if Dcondition1[i+1] == 0:
            wiewasdemol_2.append('D1!')

num_B2s = wiewasdemol_2.count('B!')
num_C2s = wiewasdemol_2.count('C!')
num_D2s = wiewasdemol_2.count('D1!')
sum_2 = num_B2s + num_C2s + num_D2s

percentage_B_2 = num_B2s / sum_2 * 100
percentage_C_2 = num_C2s / sum_2 * 100
percentage_D_2 = num_D2s / sum_2 * 100

print('F)a) Aircraft 1 number of steps condition change ', fault_1.count('Change') )
print('Aircraft 2 number of steps condition change ', fault_2.count('Change') )
print('F)b) Aircraft 1, changes, B = ', percentage_B_1,'% C = ', percentage_C_1, '% D = ', percentage_D_1, '%')
print('Aircraft 2, changes, B = ', percentage_B_2,'% C = ', percentage_C_2, '% D = ', percentage_D_2, '%')