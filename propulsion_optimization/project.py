# Code for the gas turbine system. Large or small.

# Schematic

# Necessary modules

import numpy as np
import matplotlib.pyplot as plt

### Given variables

# Temperatures

T_6 = 1800 # turbine inlet temperature
T_7 = 1200 # recuperator inlet temperature
T_32 = 276 # compressor inlet temperature 
T_10 = 363 # mixing recirculation temperature
T_amb = 288 # ambient temperature

# Pressures

P_amb = 1 # ambient pressure

# Efficiencies/Effectiveness

eta_cp = 0.9 # polytropic compressor efficiency
eta_ct = 0.92 # polytropic turbine efficiency

eff_cup = 0.92 # recuperator effectiveness
eff_cir = 0.9 # recirculator effectiveness

# Pressure drops

pi_hr = 0.97
pi_cr = 0.98

pi_ic = 0.99
pi_ec = 0.99
pi_gc = 0.99

# Gamma values

gam_h = 1.3
gam_m = 1.35
gam_c = 1.4

# Bleed

bl = 0.12

# Heat and COP

q_f = 45000000
cop = 0.7

# Recirculation ratio

phi = 0.9

##############################
# Temperature track          #
##############################

T_1 = T_amb

