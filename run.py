#########################################################################################################################
#################### ATOMICALLY RESOLVED EFFECTIVE RESISTOR NETWORK (AR-ERN) ############################################
#########################################################################################################################
#
# Version: 2.0-β
#
#########################################################################################################################


#################### PREAMBLE ###########################################################################################
#
#
import numpy as np
import math
import time
import datetime
import math
import os
#
#
#########################################################################################################################


#################### INPUTS #############################################################################################
#
#
coordinate_file = "input_example/structure_example.xyz.pe"	
						# Input file with atom coordinates and, possibly, per atom
						#     potential energy. It determines the unit for lengths.
mode='full calculation'				# To perform only part of the calculation, set:
						#    'grid' (only creates the grid for the given structure)
						#    'phase' (creates grid and assign phase to cells)
						#    'conductance' (creates grid, assign phase and conductance)
						#    'full calculation' (performed the entire routine)
grid_size=[13,18,8]				# Number of cells in each direction. Only integers allowed.
phase_assignment_method = "density"		# Cell phase can be assigned either according to the cell density
						#     (input: 'density') or the cell potential energy (input: 
						#     'potential energy').
phase_assignment_threshold = 0.04 		# Threshold for phase assignment. In the case of 'density' method, 
			   #0.039 p.e.		#     it represents the least atomic density for which the cell is 
			   #0.04 density		#     considered in pristine crystalline phase; below that 
						#     it is considered in amorphous phase.
						#     In the case of 'potential energy' method is an adimensional
						#     quantity with which the (standard deviation/mean)
						#     atomic potential energy of a cell is compared;
						#     if above that, the cell is declared in amorphous state,
						#     otherwise is crystalline.
mean_free_path={
    'a':40, 					# Electron mean free path in amorphous phase.
    'c':40					# Electron mean free path in crystal phase.
	}
#
def G_ballistic(phase,section,length):		# Formula for conductance for ballistic transport.
    if phase == 'c':
        output = (1.5*10**(-6)*section+3.87*10**(-6))*length 
    elif phase == 'a':
        output = G_ballistic('c',section,length)/2
    elif phase == 'v':
        output = 0
    else:
        print('WARNING: THERE IS A PROBLEM WITH THE PHASES')
    return output
#
def G_classical(phase,section,length):		# Formula for conductance for classical transport.
    if phase == 'c':
        output = 0.0041*section/length
    elif phase == 'a':
        output = G_classical('c',section,length)/5
    elif phase == 'v':
        output = 0
    else:
        print('WARNING: THERE IS A PROBLEM WITH THE PHASES')
    return output
#
Volt=1						# Applied Voltage. Its value does not affect in any way
						#     the calculation, apart from a scaling factor on the
						#     node voltage output.
electrode_axis = 'x'				# Axis along which electrodes are placed. Values: 'x' of 'y'.
how_to_initialize_V = 'ramp'			# Profile of V(i,j,k) at beginning of iteration procedure. Values:
						#    'ramp' for systems mostly filled with matter,
						#    'steps' for systems mostly empty.
cpu_to_use = 4					# Number of CPU to use in parallel calculation.
						#    To exploit all available cores, set very high number 
						#    and it will automatically set to max available CPUs.
iterations=1000					# Number of iterations for iterative procedure to solve 
						#    circuit equations.
#
#
#########################################################################################################################


#################### EVALUATION #########################################################################################
#
#
exec(open("src/main.py").read())
#
#
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
