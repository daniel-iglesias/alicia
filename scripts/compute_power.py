#!/usr/bin/env python

#   Copyright (C) 2007-2017 by Daniel Iglesias                            
#   Copyright (c) 2013-2018, European Atomic Energy Community (EURATOM)   
#   https://github.com/daniel-iglesias/alicia                             
#                                                                         
#   This program is free software; you can redistribute it and/or modify  
#   it under the terms of the GNU Lesser General Public License as        
#   published by the Free Software Foundation; either version 2 of the    
#   License, or (at your option) any later version.                       
#                                                                         
#   This program is distributed in the hope that it will be useful,       
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
#   GNU General Public License for more details.                          
#                                                                         
#   You should have received a copy of the GNU Lesser General Public      
#   License along with this program; if not, write to the                 
#   Free Software Foundation, Inc.,                                       
#   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             

import sys
import math
import numpy as np
from scipy.interpolate import interp1d

# Tile-dependent parameters
# read the radius from the centre of the machine
radius_in = np.loadtxt("stackA_radius.txt")

# Inputs from ALICIA
data = np.loadtxt("flux_interface.txt")
profile = data[0,1:]
time = data[1:,0]
flux = data[1:,1:]

# Inputs from FLUSH
# TODO

# Calculation of TWF
# TODO

# Calculate node_weights using profile info
node_weights = np.empty(profile.shape)
for i in range(profile.size-1) :
    node_weights[i] += (profile[i+1]-profile[i]) / 2.
    node_weights[i+1] += (profile[i+1]-profile[i]) / 2.

# Power integrated poloidally
power = np.dot(flux, node_weights)

time_power = np.column_stack(
	(time,
	power))
np.savetxt("power_local.txt",time_power)

# Total power
interp_r_fun = interp1d(radius_in[0],radius_in[1])
radius = interp_r_fun(profile)

# Calculate node_weights using radius info
node_par_weights = np.empty(profile.shape)
for i in range(radius.size-1) :
    # weight = 2*pi*r * delta_r
    node_par_weights[i] += 2*math.pi*radius[i] * (radius[i+1]-radius[i]) / 2.
    node_par_weights[i+1] += 2*math.pi*radius[i] * (radius[i+1]-radius[i]) / 2.

power_total = np.dot(flux*TWF, node_par_weights)

time_power = np.column_stack(
	(time,
	power_total))
np.savetxt("power_total.txt",time_power)

# And the energy evolution
delta_t = (time[-1]-time[0])/(time.size-1)
energy = np.zeros(power_total.shape)

for i in range(power_total.size-1) :
    # Trapezoidal rule
    energy[i+1] = energy[i] + delta_t * (power_total[i] + power_total[i+1])/2.
    
time_energy = np.column_stack(
	(time,
	energy
	))
np.savetxt("energy.txt",time_energy)

