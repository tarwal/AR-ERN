#########################################################################################################################
#################################### Conductance assignment #############################################################
#########################################################################################################################
#
print_log("**************************** Conductance assignment *****************************************************")
#


#########################################################################################################################
# BALLISTIC REGIME
#########################################################################################################################
#
# We follow the prescription explained in PHYSICAL REVIEW MATERIALS 5, 126001 (2021)
#
def run_ballistic(pacex,pacey,pacez,cellphase,numofcellsalongx,numofcellsalongy,numofcellsalongz,mean_free_path,RLandauer):
    print_log(" ")
    print_log("Ballistic transport invoked.")
    print_log(" ")
    print_log("    Grid pace(s) is smaller than electron mean free path in system.")
    print_log(" ")
    print_log("""    Resistances are being assigned following algorithm for ballistic regime... 
    (see Phys. Rev. Materials 5, 126001 (2021) for more details)""")
    print_log(" ")    

    def coresamplealongx(ny,nz):
        return [row[ny][nz] for row in cellphase]
    def coresamplealongy(nx,nz):
        return [row[nz] for row in cellphase[nx]]
    def coresamplealongz(nx,ny):
        return [row for row in cellphase[nx][ny]]
    
    samephasealongx = np.zeros((numofcellsalongx,numofcellsalongy,numofcellsalongz))
    samephasealongy = np.zeros((numofcellsalongx,numofcellsalongy,numofcellsalongz))
    samephasealongz = np.zeros((numofcellsalongx,numofcellsalongy,numofcellsalongz))

    def calcsamephasealongx(nx,ny,nz):
        sample = coresamplealongx(ny,nz)
        position = nx
        forward = 0
        while forward <= mean_free_path[cellphase[nx][ny][nz]]/pacex and position + forward + 1 <= len(sample) and sample[position+forward]==sample[position]:
            forward += 1
        if forward > 0:
            forward = forward - 1
        backward = 0
        while backward <= mean_free_path[cellphase[nx][ny][nz]]/pacex and position - backward + 1 >= 0 and sample[position-backward]==sample[position]:
            backward += 1
        if backward > 0:
            backward = backward -1
        if forward + backward + 1 < mean_free_path[sample[position]]/pacex:
            output = forward + backward + 1
        else:
            output = mean_free_path[sample[position]]/pacex 
        return output

    def calcsamephasealongy(nx,ny,nz):
        sample = coresamplealongy(nx,nz)
        position = ny
        forward = 0
        while forward <= mean_free_path[cellphase[nx][ny][nz]]/pacey and position + forward + 1 <= len(sample) and sample[position+forward]==sample[position]:
            forward += 1
        if forward > 0:
            forward = forward - 1
        backward = 0
        while backward <= mean_free_path[cellphase[nx][ny][nz]]/pacey and position - backward + 1 >= 0 and sample[position-backward]==sample[position]:
            backward += 1
        if backward > 0:
            backward = backward -1
        if forward + backward + 1 < mean_free_path[sample[position]]/pacey:
            output = forward + backward + 1
        else:
            output = mean_free_path[sample[position]]/pacey 
        return output

    def calcsamephasealongz(nx,ny,nz):
        sample = coresamplealongz(nx,ny)
        position = nz
        forward = 0
        while forward <= mean_free_path[cellphase[nx][ny][nz]]/pacez and position + forward + 1 <= len(sample) and sample[position+forward]==sample[position]:
            forward += 1
        if forward > 0:
            forward = forward - 1
        backward = 0
        while backward <= mean_free_path[cellphase[nx][ny][nz]]/pacez and position - backward + 1 >= 0 and sample[position-backward]==sample[position]:
            backward += 1
        if backward > 0:
            backward = backward -1
        if forward + backward + 1 < mean_free_path[sample[position]]/pacez: 
            output = forward + backward + 1
        else:
            output = mean_free_path[sample[position]]/pacez 
        return output

    for nx in range(numofcellsalongx):
        for ny in range(numofcellsalongy):
            for nz in range(numofcellsalongz):
                if cellphase[nx][ny][nz]=='v':
                    None
                else:
                    samephasealongx[nx,ny,nz]= calcsamephasealongx(nx,ny,nz)

    for nx in range(numofcellsalongx):
        for ny in range(numofcellsalongy):
            for nz in range(numofcellsalongz):
                if cellphase[nx][ny][nz]=='v':
                    None
                else:
                    samephasealongy[nx,ny,nz]= calcsamephasealongy(nx,ny,nz)

    for nx in range(numofcellsalongx):
        for ny in range(numofcellsalongy):
            for nz in range(numofcellsalongz):
                if cellphase[nx][ny][nz]=='v':
                    None
                else:
                    samephasealongz[nx,ny,nz]= calcsamephasealongz(nx,ny,nz)

    for nx in range(numofcellsalongx):
        for ny in range(numofcellsalongy):
            for nz in range(numofcellsalongz):
                if cellphase[nx][ny][nz] == 'v':
                    gx[nx,ny,nz]=0 
                else:
                    gx[nx,ny,nz]=G_ballistic(cellphase[nx][ny][nz],pacey*pacez,samephasealongx[nx, ny, nz])

    for nx in range(numofcellsalongx):
        for ny in range(numofcellsalongy):
            for nz in range(numofcellsalongz):
                if cellphase[nx][ny][nz] == 'v':
                    gy[nx,ny,nz]=0 
                else:
                    gy[nx,ny,nz]=G_ballistic(cellphase[nx][ny][nz],pacex*pacey,samephasealongy[nx, ny, nz])

    for nx in range(numofcellsalongx):
        for ny in range(numofcellsalongy):
            for nz in range(numofcellsalongz):
                if cellphase[nx][ny][nz] == 'v':
                    gz[nx,ny,nz]=0 
                else:
                    gz[nx,ny,nz]=G_ballistic(cellphase[nx][ny][nz],pacex*pacey,samephasealongz[nx, ny, nz])

    print_log("    All cells have now associated a resistance value.")
    print_log(" ")
    print_log("Effective Resistor Network created.")
    print_log(" ")
    return gx,gy,gz
#
#
#
#########################################################################################################################
# CLASSICAL REGIME
#########################################################################################################################
#
#
def run_classical(pacex,pacey,pacez,cellphase,numofcellsalongx,numofcellsalongy,numofcellsalongz,mean_free_path,RLandauer):
    print_log(" ")
    print_log("Classical Effective Resistor Network.")
    print_log(" ")
    print_log("    Grid pace(s) is larger than electron mean free path in system.")
    print_log(" ")
    print_log("""    Cell conductances are being assigned...""")
    print_log(" ")    

    for nx in range(numofcellsalongx):
        for ny in range(numofcellsalongy):
            for nz in range(numofcellsalongz):
                if cellphase[nx][ny][nz] == 'v':
                    gx[nx,ny,nz]=0 
                    gy[nx,ny,nz]=0 
                    gz[nx,ny,nz]=0 
                else:
                    gx[nx,ny,nz]=G_classical(cellphase[nx][ny][nz],pacez*pacey,pacex)
                    gy[nx,ny,nz]=G_classical(cellphase[nx][ny][nz],pacez*pacex,pacey)
                    gz[nx,ny,nz]=G_classical(cellphase[nx][ny][nz],pacex*pacey,pacez)
    
    print_log("    All cells have now associated a resistance value.")
    print_log(" ")
    print_log("Effective Resistor Network created.")
    print_log(" ")
    return gx,gy,gz
#
#
#########################################################################################################################


#########################################################################################################################
# ASSIGN CONDUCTANCES
#########################################################################################################################
#
#
tmp_list=list(mean_free_path.values())
tmp_list.sort()
least_mean_free_path=tmp_list[0]
del tmp_list
#
# Initialize cell condunctances
gx = np.zeros((numofcellsalongx,numofcellsalongy,numofcellsalongz))
gy = np.zeros((numofcellsalongx,numofcellsalongy,numofcellsalongz))
gz = np.zeros((numofcellsalongx,numofcellsalongy,numofcellsalongz))
#
# ASSIGNING RESISTANCE TO CELLS
#
if pacex<least_mean_free_path or pacey<least_mean_free_path or pacez<least_mean_free_path:
    gx,gy,gz = run_ballistic(pacex,pacey,pacez,cellphase,numofcellsalongx,numofcellsalongy,numofcellsalongz,mean_free_path,G_ballistic)
else:
    gx,gy,gz = run_classical(pacex,pacey,pacez,cellphase,numofcellsalongx,numofcellsalongy,numofcellsalongz,mean_free_path,G_ballistic)
#
#
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

