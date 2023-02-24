import numpy as np
import matplotlib.pyplot as plt

# Provided
box_length = 1000 # meters
box_volume = box_length**3 # cubic meters

Navo = 6.02e26 # molec/kmol
rho = 1.0 # kg/m^3
Mair = 29.0 # kg/kmol
P_ocean = 0.0 # kgN/day
P_land = 0.1 # kgN/day
VMR = 10/1e12
timestep = 3600 # seconds so one hour
L = 1/216000 # 1/s
T = 293 # K

VMR_CO = 140/1e9 # ppbv
CO_concentration_initialvalue = [VMR_CO * Navo * rho / (Mair * 1e6)] # molec/cm^3

VMR_CH4 = 1900/1e9 # ppbv
CH4_concentration_initialvalue = [VMR_CO * Navo * rho / (Mair * 1e6)] # molec/cm^3

cair = rho / Mair * Navo * 10e3 * (1/10e6) # molec/m^3

k4 = 1.75 * 1.57e-13 + cair * 3.54e-33
k5 = 1.85e-20 * np.exp(2.82*np.log(T)-987/T)
k7 = 1e-14 * np.exp(-490/T)
k8 = 1.7e-12 * np.exp(-940/T)

#CO production
P4 = k4 * CO_concentration * OH_concentration

#CH4 production
P5 = k5 * CH4_concentration * OH_concentration

#O3 production
P7 = k7 * O3_concentration * HO2_concentration
P8 = k8 * O3_concentration * OH_concentration