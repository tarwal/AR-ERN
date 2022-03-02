#########################################################################################################################
#################### Coordinate file importing and creation of the grid #####################################################
#########################################################################################################################
#
# xyz file reading adapted from:
#    
print_log("**************************** File import and grid creation *********************************************")
print_log(" ")
#
#  Importing file .xyz
#
if coordinate_file[-4::] == '.xyz':
    print_log("Reading file " + coordinate_file + " in XYZ format...")
    print_log(" ")
    print_log("    Check on the number of coordinates and the number of atoms declared in the coordinate file.")
    print_log(" ")
    #
    atomtypes = []
    coordinates = []
    xyz = open(coordinate_file)
    n_atoms = int(xyz.readline())
    xyz.readline()
    relevant_lines = range(n_atoms)
    for pos, line in enumerate(xyz):
        if pos in relevant_lines:      
            atom,x,y,z = line.split()  # Be careful with the file formatting: 
				     # in this line I specify that the elements of the same row are 
				     # separated by a single space. You can change this 
				     # by changing the argument of line.split().
            atomtypes.append(atom)
            coordinates.append([float(x), float(y), float(z)])
    xyz.close()
    if n_atoms != len(coordinates):
        raise ValueError("File says %d atoms but read %d points." % (n_atoms, len(coordinates)))
    #
    print_log("Atomic coordinates correctly extracted from coordinate file.")
    print_log(" ")
else:
    None
#
#  Alternatively, importing file .xyz.pe.dat
#
if coordinate_file[-7::] == '.xyz.pe':
    print_log("Reading file " + coordinate_file + " in XYZ+PE format...")
    print_log(" ")
    print_log("    Check on the number of coordinates and the number of atoms declared in the coordinate file.")
    print_log(" ")
    #
    atomtypes = []
    coordinates = []
    pot_ene = []
    linenumber = 1
    xyzpe = open(coordinate_file)##+'.dat')
    xyzpe.readline()
    xyzpe.readline()
    xyzpe.readline()
    n_atoms = int(xyzpe.readline())
    xyzpe.readline()
    xyzpe.readline()
    xyzpe.readline()
    xyzpe.readline()
    xyzpe.readline()
    relevant_lines = range(n_atoms)    
    for pos, line in enumerate(xyzpe):
        if pos in relevant_lines:      
            atom,x,y,z,pe = line.split()  # Be careful with the file formatting: 
				   # in this line I specify that the elements of the same row are 
				   # separated by a single space. You can change this 
				   # by changing the argument of line.split().
            atomtypes.append(atom)
            pot_ene.append(float(pe))
            coordinates.append([float(x), float(y), float(z)])
    xyzpe.close()
    if n_atoms != len(coordinates):
        raise ValueError("File says %d atoms but read %d points." % (n_atoms, len(coordinates)))
    #
    print_log("Atomic coordinates correctly extracted from coordinate file.")
    print_log(" ")
else:
    None
#
#  If the file is not in one of the above formats, abort evaluation.
#
if 'coordinates' in locals():
    None
else:
    import sys
    sys.exit('''ERROR MESSAGE: COORDINATE FILE NOT IMPORTED. 
    PLEASE MAKE SURE THE FILE IS IN THE RIGHT FOLDER AND HAS THE RIGHT FORMATTING.
    EVALUATION ABORTED''')
#
#
#  Grid creation.
#
print_log("Finding system corners...")
#
numofatoms=len(coordinates)
atomcoordinates=np.array(coordinates)
xminimo = min(atomcoordinates[:,0])
yminimo = min(atomcoordinates[:,1])
zminimo = min(atomcoordinates[:,2])
xmassimo = max(atomcoordinates[:,0])
ymassimo = max(atomcoordinates[:,1])
zmassimo = max(atomcoordinates[:,2])
#
print_log(" ")
print_log("    Box opposite corners:")
print_log("        ",[xminimo,yminimo,zminimo])
print_log("        ",[xmassimo,ymassimo,zmassimo])
print_log(" ")
print_log("Dividing the space in cells...")
#
pacex = (xmassimo - xminimo)/grid_size[0]
pacey = (ymassimo - yminimo)/grid_size[1]
pacez = (zmassimo - zminimo)/grid_size[2]
numofcellsalongx = grid_size[0]
numofcellsalongy = grid_size[1]
numofcellsalongz = grid_size[2]
#
print_log(" ")
print_log("    Grid size: ")
print_log("        ", numofcellsalongx,"x", numofcellsalongy,"x", numofcellsalongz)
print_log("    Cell size: ")
print_log("        ",pacex,"x", pacey,"x", pacez)
print_log(" ")
print_log("Assigning indexing to cells...")
#
def calculatecelllimits(nx,ny,nz):
    return [[xminimo+nx*pacex,xminimo+(nx+1)*pacex],
 		[yminimo+ny*pacey,yminimo+(ny+1)*pacey],
		[zminimo+nz*pacez,zminimo+(nz+1)*pacez]]
cellcoordinates=np.array([[[
	calculatecelllimits(nx,ny,nz)
	for nz in range(0,numofcellsalongz)]
	for ny in range(0,numofcellsalongy)]
	for nx in range(0,numofcellsalongx)]
	)
#
print_log(" ")
print_log("""    First cell coordinates (format: [[xmin,xmax],[ymin,ymax],[zmin,zmax]]):
    """+str(cellcoordinates[0,0,0]).replace('\n','\n\t'))
print_log("""    Last cell coordinates:
    """+str(cellcoordinates[numofcellsalongx-1, numofcellsalongy-1, numofcellsalongz-1]).replace('\n','\n\t'))
print_log(" ")
#
#
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################