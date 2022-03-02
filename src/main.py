

start_time = time.time()


#**************************** Log files initialization **************************************************
#
# 
i=1
while coordinate_file[-i] != '/' and i<len(coordinate_file):
	i+=1
if coordinate_file[-i] == '/':
    if coordinate_file[-4::] == '.xyz':
        coordinate_file_no_ext=coordinate_file[-i+1:].replace('.xyz','_')
    elif coordinate_file[-7::] == '.xyz.pe':
        coordinate_file_no_ext=coordinate_file[-i+1:].replace('.xyz.pe','_')
    else:
        print('WARNING: THERE IS A PROBLEM WITH THE EXTENSION OF THE STRUCTURE FILE.')
else:
    if coordinate_file[-4::] == '.xyz':
        coordinate_file_no_ext=coordinate_file[-i:].replace('.xyz','_')
    elif coordinate_file[-7::] == '.xyz.pe':
        coordinate_file_no_ext=coordinate_file[-i:].replace('.xyz.pe','_')
    else:
        print('WARNING: THERE IS A PROBLEM WITH THE EXTENSION OF THE STRUCTURE FILE.')
del i
#
log_file_name =coordinate_file_no_ext+'AR-ERN.log'
if os.path.exists(log_file_name):
    os.remove(log_file_name)
log_file = open(log_file_name, "a")
log_file.write("""  AR-ERN INPUT FILE

""")
log_file.close()
f1 = log_file_name
f2 = "run.py"
f1 = open(f1, 'a+')
f2 = open(f2, 'r')
f1.write(f2.read())
f1.write("""

  AR-ERN LOG FILE

""")
f1.close()
f2.close()
#
#
def convertTuple(tup):
    if type(tup) == str:
        return tup
    else:
        st = ''.join(map(str, tup))
        return st
def print_log(*x):
    print(convertTuple(x))
    log_file = open(log_file_name, "a")
    log_file.write('\n')
    log_file.write(convertTuple(x))
    log_file.close()
#
#
print_log("""
========================================================================================================
=================================== AR-ERN INITIATED... ================================================
========================================================================================================
""")
#
#
#********************************************************************************************************



#**************************** Output files initialization ***********************************************
#
# 
if os.path.isdir('AR-ERN_output'):
    None
else:        
    os.mkdir('AR-ERN_output')
#
output_file_name_Resistance='AR-ERN_output/'+coordinate_file_no_ext+'Resistance.csv'
output_file_name_phases='AR-ERN_output/'+coordinate_file_no_ext+'phases.csv'
output_file_name_conductances='AR-ERN_output/'+coordinate_file_no_ext+'conductances.csv'
output_file_name_voltages='AR-ERN_output/'+coordinate_file_no_ext+'voltages.csv'
#
if os.path.exists(output_file_name_Resistance) or os.path.exists(output_file_name_phases) or os.path.exists(output_file_name_conductances) or os.path.exists(output_file_name_voltages):
    if os.path.exists(output_file_name_Resistance):
        os.remove(output_file_name_Resistance)
    if os.path.exists(output_file_name_phases):
        os.remove(output_file_name_phases)
    if os.path.exists(output_file_name_conductances):
        os.remove(output_file_name_conductances)
    if os.path.exists(output_file_name_voltages):
        os.remove(output_file_name_voltages)
    print_log('    Output files with same name already existing. Files replaced with new ones.')
    print_log(' ')
else:
    print_log('    New output files created.')
    print_log(' ')

#
output_file_Resistance = open(output_file_name_Resistance, "a")
output_file_Resistance.write("Iteration,Resistance")
output_file_Resistance.close()
output_file_phases = open(output_file_name_phases, "a")
output_file_phases.write("Cell index i,Cell index j,Cell index k,Cell phase")
output_file_phases.close()
output_file_conductances = open(output_file_name_conductances, "a")
output_file_conductances.write("Cell index i,Cell index j,Cell index k,Cell conductance along x,Cell conductance along y,Cell conductance along z")
output_file_conductances.close()
output_file_voltages = open(output_file_name_voltages, "a")
output_file_voltages.write("Node index i,Node index j,Node index k,Node voltage")
output_file_voltages.close()
#
#
#********************************************************************************************************


#**************************** File import and grid creation *********************************************
#
# Requires: coordinate_file, grid_size
# Returns: atomtypes, atomcoordinates, cellcoordinates, numofatoms, numofcellsalongx, numofcellsalongy, 
#          numofcellsalongz, pacex, pacey, pacez, (pot_ene)
#
exec(open("src/create_grid.py").read())
#
#********************************************************************************************************


#**************************** Phase assignment **********************************************************
#
# Requires: pacex, pacey, pacez, phase_assignment_threshold, cellcoordinates, atomcoordinates, 
#           numofcellsalongx, numofcellsalongy, numofcellsalongz, numofatoms, (pot_ene)
# Returns: cellphase
#
if mode == 'grid':
    None
else:
    if phase_assignment_method == "density":
        exec(open("src/assign_phase__density.py").read())
    else:
        exec(open("src/assign_phase__pot_ene.py").read())
    outputfile = open(output_file_name_phases, 'a')
    for i in range(numofcellsalongx):
        for j in range(numofcellsalongy):
            for k in range(numofcellsalongz):
                outputfile.write("\n"+str(i)+","+str(j)+","+str(k)+","+str(cellphase[i][j][k]))
    outputfile.close()
#
#********************************************************************************************************


#**************************** Conduction assignment *****************************************************
#
# Requires: pacex, pacey, pacez, cellphase, numofcellsalongx, numofcellsalongy, numofcellsalongz,
#           mean_free_path, G_ballistic
# Returns: gx,gy,gz
#
if mode == 'grid' or mode == 'phase':
    None
else:
    exec(open("src/assign_conductance.py").read())
    outputfile = open(output_file_name_conductances, 'a')
    for i in range(numofcellsalongx):
        for j in range(numofcellsalongy):
            for k in range(numofcellsalongz):
                outputfile.write("\n"+str(i)+","+str(j)+","+str(k)+","+str(gx[i,j,k])+","+str(gy[i,j,k])+","+str(gz[i,j,k]))
    outputfile.close()
#
#********************************************************************************************************


#**************************** Solve circuit equations ***************************************************
#
# Requires: coordinate_file, cpu_to_use, how_to_initialize_V, electrode_axis, gx, gy, gz, Volt, iterations,
#           cellphase, numofcellsalongx, numofcellsalongy, numofcellsalongz,
# Returns: v, Rlist
#
if mode == 'full calculation':
    exec(open("src/solve_circuit_equations.py").read())
    outputfile = open(output_file_name_voltages, 'a')
    for i in range(numofcellsalongx):
        for j in range(numofcellsalongy):
            for k in range(numofcellsalongz):
                outputfile.write("\n"+str(i)+","+str(j)+","+str(k)+","+str(v(i,j,k)))
    outputfile.close()
else:
    None
#
#********************************************************************************************************

#**************************** Wrap up *******************************************************************
#
print_log(" ")
print_log('Total time: '+str(datetime.timedelta(seconds=round(time.time() - start_time))))
print_log("""
========================================================================================================	
=================================== ...AND CLOSED ======================================================
========================================================================================================
""")
#
#
#########################################################################################################################
