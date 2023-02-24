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

# Could be wrong
L = 1/216000 # 1/s

# Timestep
time = np.arange(0, 30*24*3600, timestep) # 30 days
NOx_concentration = [10 * VMR * Navo * rho / (Mair * 1e6)] # molec/cm^3

# Initialize NOx_total array with zeros
NOx_total = np.zeros(len(time))
NOx_total[0] = NOx_concentration[0] * box_volume

# Numerical solution using semi-implicit Euler scheme
for i in range(len(time)-1):
    P = P_ocean if i*timestep <= 10*24*3600 else P_land
    NOx_total[i+1] = (NOx_total[i] + P*timestep) / (1 + L*timestep)
    NOx_concentration.append(NOx_total[i+1] / box_volume)

# Plot results
plt.plot(time/(24*3600), np.array(NOx_concentration)*box_volume)
plt.xlabel('Time (days)')
plt.ylabel('NOx (kgN)')
plt.title('Temporal Evolution of NOx over 30 days')
plt.show()