
#########################################################################################################################
############################# Phase Assignment ##########################################################################
#########################################################################################################################
#
#
print_log("**************************** Phase assignment **********************************************************")
#
# ASSIGNING PHASE TO GRID CELLS
#
# Method based on atomic potential energy.
#
if 'pot_ene' in locals():
    None
else:
    import sys
    sys.exit('''ERROR MESSAGE: ATOMIC POTENTIAL ENERGY NOT IMPORTED. 
    PLEASE MAKE SURE THAT THE COORDINATE FILE CONTAINS IT.
    EVALUATION ABORTED''')
print_log(" ")
print_log("Computing mean and standard deviation of atomic poential energy for each cell...")
print_log(" ")
#
cellxmin = cellcoordinates[0,0,0][0,0]
cellymin = cellcoordinates[0,0,0][1,0]
cellzmin = cellcoordinates[0,0,0][2,0]
atomcoordinates=np.array(atomcoordinates)
poteneincell = [[[[] for _ in range(numofcellsalongz)] for _ in range(numofcellsalongy)] for _ in range(numofcellsalongx)]
atomsincell = np.zeros((numofcellsalongx, numofcellsalongy, numofcellsalongz))

for i in range(numofatoms):
    position_x=math.floor((atomcoordinates[i,0]-cellxmin)/pacex)
    position_y=math.floor((atomcoordinates[i,1]-cellymin)/pacey)
    position_z=math.floor((atomcoordinates[i,2]-cellzmin)/pacez)
    if position_x == numofcellsalongx or position_y == numofcellsalongy or position_z == numofcellsalongz:
        None
    else:
        atomsincell[position_x,position_y,position_z] += 1 
        poteneincell[position_x][position_y][position_z].append(pot_ene[i])
if numofatoms==np.sum(atomsincell):
    test_result = 'Yes.'
else:
    test_result = '    No, '+str(int(numofatoms-np.sum(atomsincell)))+''' atoms where on the farthest border of the grid and have been neglected.
    If such a number looks to you too high in comparison with the total number of particles ('''+str(numofatoms)+'''), 
    then modify the way the grid is defined.'''
print_log("    Have all atoms been assigned to a cell? ")
print_log(" ")
print_log(test_result)
print_log(" ")
print_log("Assigning phase based on atomic potential energy...")
print_log(" ")
#
cellphase = np.zeros((numofcellsalongx, numofcellsalongy, numofcellsalongz))
cellphase = cellphase.tolist()
for nx in range(numofcellsalongx):
    for ny in range(numofcellsalongy):
        for nz in range(numofcellsalongz):
            if atomsincell[nx][ny][nz]==0:
                cellphase[nx][ny][nz]='v'
            elif np.absolute(np.std(np.array(poteneincell[nx][ny][nz]))/np.mean(np.array(poteneincell[nx][ny][nz]))) < phase_assignment_threshold:
                cellphase[nx][ny][nz]='c'
            else:
                cellphase[nx][ny][nz]='a'
#
cellphaseflatten = sum(cellphase, [])
cellphaseflatten = sum(cellphaseflatten, [])
if len(cellphaseflatten)==(numofcellsalongx*numofcellsalongy*numofcellsalongz):
    test_result = 'Yes.'
else:
    test_result = '******************* No! ERROR, ABORT MANUALLY, PLEASE! *******************'
#
print_log("    Does the number of cells with assigned phase equals the total number of cells? ",test_result)
print_log(" ")
#
vacuumpercentage=round(100*(sum(x.count('v') for x in cellphaseflatten)/len(cellphaseflatten)),1)
crystalpercentage=round(100*(sum(x.count('c') for x in cellphaseflatten)/len(cellphaseflatten)),1)
amorphouspercentage=round(100*(sum(x.count('a') for x in cellphaseflatten)/len(cellphaseflatten)),1)
#
print_log("    Phase percentages: ", amorphouspercentage,"% of amorphous phase, ", crystalpercentage,"% of crystal phase, ", vacuumpercentage,"% of vacuum.")
print_log(" ")
print_log("Phases assigned.")
print_log(" ")
#
#
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
