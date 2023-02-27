import numpy as np
import matplotlib.pyplot as plt

# Provided
box_length = 1000 # meters
box_volume = box_length**3 # cubic meters
box_to_cm = 1e15
Navo = 6.02e26 # molec/kmol
rho = 1.0 # kg/m^3
Mair = 29.0 # kg/kmol
P_ocean = 0.0 # kgN/day
P_land = 0.1 # kgN/day
VMR = 10/1e12 # pptv
timestep = 3600 # seconds so one hour
T = 293 # K
M_N = 14.0 # kg/kmol
M_H = 1 # kg/kmol
M_C = 12 # kg/kmol
M_O = 16 # kg/kmol
P_CO_land = 12
P_CH4_land = 0.24

L = 1/216000 # 1/s

# Timestep
time = np.arange(0, 30*24*3600, timestep) # 30 days

cair = rho / Mair * Navo * (1/1e6) # molec/cm^3

k1 = 3e-12 * np.exp(-1500/T)
k2 = 5e-3
k3 = 5.1e-12 * np.exp(210/T)
k4 = 1.57e-13 + cair * 3.54e-33
k5 = 1.85e-20 * np.exp(2.82*np.log(T)-987/T)
k6 = 3.3e-12 * np.exp(270/T)
k7 = 1e-14 * np.exp(-490/T)
k8 = 1.7e-12 * np.exp(-940/T)

# CO intial
VMR_CO = 140/1e9 # ppbv
CO_concentration_0 = VMR_CO * Navo * rho / (Mair * 1e6) # molec/cm^3


# CH4 initial
VMR_CH4 = 1900/1e9 # ppbv
CH4_concentration_0 = VMR_CH4 * Navo * rho / (Mair * 1e6) # molec/cm^3

# O3 initial
VMR_O3 = 45/1e9 # ppbv
O3_concentration_0 = VMR_O3 * Navo * rho / (Mair * 1e6) # molec/cm^3

# OH constant
OH_concentration_constant = 1.5e6 # molec/cm^3

# HO2 constant
HO2_concentration_constant = 450e6 # molec/cm^3

# O initial
O_concentration_0 = O3_concentration_0 * 5e-6

# CO loss
# Initialize CO array with zeros
CO_concentration = np.zeros(len(time))
CO_concentration[0] = CO_concentration_0
# Initialize CO array with zeros ppbv

#Production term
P4 = P_CO_land * Navo / ((M_C+M_O) * 24 * timestep * box_to_cm)


# Numerical solution using semi-implicit Euler scheme
for i in range(len(time)-1):
    L4 = k4 * OH_concentration_constant
    if i*timestep <= 10*24*3600:
        P_CO = 0
    else:
        P_CO = P4
    CO_concentration[i+1] = (CO_concentration[i] + P_CO * timestep) / (1 + L4*timestep)

# CH4 loss
# Initialize CH4 array with zeros
CH4_concentration = np.zeros(len(time))
CH4_concentration[0] = CH4_concentration_0

#Production term
P5 = P_CH4_land * Navo / ((M_C+(M_H*4)) * 24 * timestep * box_to_cm)

# Numerical solution using semi-implicit Euler scheme
for i in range(len(time)-1):
    L5 = k5 * OH_concentration_constant
    if i*timestep <= 10*24*3600:
        P_CH4 = 0
    else:
        P_CH4 = P5
    CH4_concentration[i+1] = (CH4_concentration[i] + P_CH4 * timestep) / (1 + L5 *timestep)

# Hier begint de herres: NOx and O3 production and loss

P_land_new = P_land * Navo / (M_N * 24 * timestep * box_to_cm)

NOx_concentration_0 = VMR * Navo * rho / (Mair * 1e6) # molec/cm^3

# Initialize NOx_total array with zeros
NOx_concentration = np.zeros(len(time))
NOx_concentration[0] = NOx_concentration_0

# Initialize O3, O and NO array with zeros
O3_concentration = np.zeros(len(time))
O3_concentration[0] = O3_concentration_0

O_concentration = np.zeros(len(time))
O_concentration[0] = O_concentration_0

NO_concentration_0 = NOx_concentration_0 * ((k2 + k3 * O_concentration_0)/(k2 + k3 * O_concentration_0 + k1 * O3_concentration_0))
NO_concentration = np.zeros(len(time))
NO_concentration[0] = NO_concentration_0

NO2_concentration_0 = NO_concentration_0 * O3_concentration_0 * k1 / (k2+k3 * O_concentration_0)
NO2_concentration = np.zeros(len(time))
NO2_concentration[0] = NO2_concentration_0
L9_0 = k1 * NO_concentration_0
L9 = np.zeros(len(time))
L9[0] = L9_0


P_O3 = np.zeros(len(time))
P_O3[0] = k6 * HO2_concentration_constant * NO_concentration[0]

L7 = k7 * HO2_concentration_constant
L8 = k8 * OH_concentration_constant
LO3 = L7+L8

# Numerical solution using semi-implicit Euler scheme
for i in range(len(time)-1):
    if i*timestep <= 10*24*3600:
        P = P_ocean
    else:
        P = P_land_new
    NOx_concentration[i+1] = (NOx_concentration[i] + P*timestep) / (1 + L*timestep)
    O_concentration[i+1] = O3_concentration[i] * 5e-6
    NO_concentration[i+1] = NOx_concentration[i] *(k2 + k3 * O_concentration[i])/(k2 + k3*O_concentration[i] + k1 * O3_concentration[i])
    P_O3[i+1] = k6 * HO2_concentration_constant * NO_concentration[i]
    L9[i+1] = k1 * NO_concentration[i]
    O3_concentration[i+1] = (O3_concentration[i] + P_O3[i] * timestep) / (1 + (LO3+L9[i]) *timestep)
    NO2_concentration[i+1] = NO_concentration[i] * O3_concentration[i] * k1 / (k2+k3 *O_concentration[i])

plt.plot(time/(24*3600), np.array(O3_concentration * 1e10 * rho / (1/Mair *1/1e6 *Navo)), label='O3 (0.1 ppbv)')
plt.plot(time/(24*3600), np.array(CH4_concentration* 1e8 * rho / (1/Mair *1/1e6 *Navo)), label='CH4 (10 ppbv)')
plt.plot(time/(24*3600), np.array(CO_concentration* 1e9 * rho / (1/Mair *1/1e6 *Navo)), label='CO (ppbv)')
plt.plot(time/(24*3600), np.array(NOx_concentration* 1e12 * rho / (1/Mair *1/1e6 *Navo)), label='NOx (pptv)')
plt.xlabel('Time (days)')
plt.ylabel('Concentration')
plt.title('Temporal Evolution over 30 days')
plt.legend()
plt.show()
