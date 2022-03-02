#########################################################################################################################
#################################### Solving the circuit equations ######################################################
#########################################################################################################################
#
print_log("**************************** Resolution of circuit equations *******************************************")
#
# Preamble
#
import multiprocessing as mp
import sys
import os
#
# Fire up Message
#
print_log(" ")
print_log("Solving circuit equations for Effective Resistor Network...")
#
# Preparing for parallel processing 1: setting the number of processors to use
#
maxcpu=mp.cpu_count()
if cpu_to_use > maxcpu:
    cpunumber = maxcpu
else:
    cpunumber = cpu_to_use
print_log(" ")
print_log("    Number of processors that will be used: "+str(cpunumber))
print_log(" ")
print_log("    Getting ready to parallel computation...")
#
# Preparing for parallel processing 2: index unwrapping functions
# Some more details.
# So far V(i,j,k) was a tensor. To make parallel computing work, it must be converted to a list.
# Moreover, to apply the chosen algorithm to solve the equations,
# we have to distinguish two lists to update in sequence, one with i+j+k=even and one with i+j+k=odd.
# Here I define two maps (in the form of numpy arrays), one that takes aas input i,j,k and gives the corresponding
# index of the even/odd index list, and the invers one.
#
zerobin=0
onebin=0
oddoreven=[]
for i in range(numofcellsalongx+1):
    for j in range(numofcellsalongy+1):
        for k in range(numofcellsalongz+1):
            oddoreven.append((i+j+k)%2)
fromijktoevenoddindex=[[[0 for k in range(numofcellsalongz+1)] for j in range(numofcellsalongy+1)] for i in range(numofcellsalongx+1)]
for i in range(numofcellsalongx+1):
    for j in range(numofcellsalongy+1):
        for k in range(numofcellsalongz+1):
            fulllistindex=np.ravel_multi_index((i,j,k),(numofcellsalongx+1,numofcellsalongy+1,numofcellsalongz+1))+1
            upto=oddoreven[:fulllistindex]
            num=i+j+k
            if (num % 2) ==0:
                fromijktoevenoddindex[i][j][k]=zerobin
                zerobin +=1
            else:
                fromijktoevenoddindex[i][j][k]=onebin
                onebin +=1
#
# Preparing for parallel processing 3: converting V(i,j,k) from tensor to two lists: Vodd and Veven, where the parity refers to i+j+k
#
dimV=(numofcellsalongx+1)*(numofcellsalongy+1)*(numofcellsalongz+1)
dimVodd=dimV//2
if (dimV % 2) ==0:
    dimVeven=dimVodd
else:
    dimVeven=dimVodd+1
Vodd=np.zeros(dimVodd)
Veven=np.zeros(dimVeven)
#
# Preparing for parallel processing 4: index wrapping function
#
fromevenoddindextoijk=[[0 for index in range(dimV)] for parity in range(2)]
evenbin=0
oddbin=0
for index in range(len(oddoreven)):
    if oddoreven[index]==0:
        fromevenoddindextoijk[0][evenbin]=np.unravel_index(index,(numofcellsalongx+1,numofcellsalongy+1,numofcellsalongz+1))
        evenbin +=1        
    else:
        fromevenoddindextoijk[1][oddbin]=np.unravel_index(index,(numofcellsalongx+1,numofcellsalongy+1,numofcellsalongz+1))
        oddbin +=1        
#
# Preparing for parallel processing 5: defining a function v(i,j,k) that recalls Vodd/Veven to use in the equation
#
def v(i,j,k):
    num=i+j+k
    if (num % 2) ==0:
        return Veven[fromijktoevenoddindex[i][j][k]]
    else:
        return Vodd[fromijktoevenoddindex[i][j][k]]
#
# RHS of the equation to solve (the LHS being V(i,j,k), see Eq. (1) of Phys. Rev. Materials 5, 126001 (2021))
# At i=0 or j=0 or k=0 the RHS is not defined (it would require Rx(-1,j,k) etc.) but the function I am defining 
# takes also values there (it simply takes literally the value '-1' as last element of the list gx[-1,j,k] etc.,
# in a sort of periodic boundary continuation, but in fact it will not be used with those indices.
#
den = np.zeros((numofcellsalongx,numofcellsalongy,numofcellsalongz))
for i in range(numofcellsalongx):
    for j in range(numofcellsalongy):
        for k in range(numofcellsalongz):
            den[i,j,k]=gx[i-1,j,k]+gx[i,j,k]+gy[i,j-1,k]+gy[i,j,k]+gz[i,j,k-1]+gz[i,j,k]
def rhs(i,j,k):
    tmp_den = den[i,j,k]
    if tmp_den == 0: 	    # The denominator of the RHS of the equation is 0 only when also the numerator is.
			    # In that case the entire equation is not well defined. 
	 		    # Instead of creating an exeption for those indices, we define the equation also
			    # in that case by assigning an arbitrary value to the RHS (here chosen to be 0).
			    # The result should not (and, as far as I tested) does not depend on such a value.
        return 0
    else:
        return ((1/tmp_den)*(gx[i-1,j,k]*v(i-1,j,k)+gx[i,j,k]*v(i+1,j,k)+gy[i,j-1,k]*v(i,j-1,k)+gy[i,j,k]*v(i,j+1,k)+gz[i,j,k-1]*v(i,j,k-1)+gz[i,j,k]*v(i,j,k+1)))
#
# Deciding where to put the electrodes
# Note that the system resistance is here defined as the electrodes Voltage / the net current flowing to the system (incoming - outcoming)
#
if electrode_axis == 'x':
    def calculate_R():
        return Volt/((sum((v(0,j,k)-v(1,j,k))*gx[0,j,k]+(v(numofcellsalongx-1,j,k)-v(numofcellsalongx,j,k))*gx[numofcellsalongx-1,j,k] for j in range(numofcellsalongy) for k in range(numofcellsalongz)))/2)
    def set_electrodes(x):
        for j in range(numofcellsalongy+1):
            for k in range(numofcellsalongz+1):
                # Positive electrode
                evenoroddindex=fromijktoevenoddindex[0][j][k]
                num=0+j+k
                if (num % 2) == 0:
                    Veven[evenoroddindex]=x/2
                else:
                    Vodd[evenoroddindex]=x/2
                # Negative electrode
                evenoroddindex=fromijktoevenoddindex[numofcellsalongx][j][k]
                num=numofcellsalongx+j+k
                if (num % 2) == 0:
                    Veven[evenoroddindex]=-x/2
                else:
                    Vodd[evenoroddindex]=-x/2
    def initialize_V(x):
        for j in range(numofcellsalongy+1):
            for k in range(numofcellsalongz+1):
                for i in range(numofcellsalongx+1):
                    evenoroddindex=fromijktoevenoddindex[i][j][k]
                    num=i+j+k
                    if (num % 2) == 0:
                        Veven[evenoroddindex]=(x/2)-(i/numofcellsalongx)*x
                    else:
                        Vodd[evenoroddindex]=(x/2)-(i/numofcellsalongx)*x
else:
    def calculate_R():
        return Volt/((sum((v(i,numofcellsalongy-1,k)-v(i,numofcellsalongy,k))*gy[i,numofcellsalongy-1,k]+(v(i,0,k)-v(i,1,k))*gy[i,0,k] for i in range(numofcellsalongx) for k in range(numofcellsalongz)))/2)
    def set_electrodes(x):
        for i in range(numofcellsalongx+1):
            for k in range(numofcellsalongz+1):
                # Positive electrode
                evenoroddindex=fromijktoevenoddindex[i][0][k]
                num=i+0+k
                if (num % 2) == 0:
                    Veven[evenoroddindex]=x/2
                else:
                    Vodd[evenoroddindex]=x/2
                # Negative electrode
                evenoroddindex=fromijktoevenoddindex[i][numofcellsalongy][k]
                num=i+numofcellsalongy+k
                if (num % 2) == 0:
                    Veven[evenoroddindex]=-x/2
                else:
                    Vodd[evenoroddindex]=-x/2
    def initialize_V(x):
        for i in range(numofcellsalongx+1):
            for k in range(numofcellsalongz+1):
                for j in range(numofcellsalongy+1):
                    evenoroddindex=fromijktoevenoddindex[i][j][k]
                    num=i+j+k
                    if (num % 2) == 0:
                        Veven[evenoroddindex]=(x/2)-(j/numofcellsalongx)*x
                    else:
                        Vodd[evenoroddindex]=(x/2)-(j/numofcellsalongx)*x

#
# Initialize Veven and Vodd with the electrodes.
# If set to 'discontinuous' the tensor V(i,j,k) is zero everywhere
# bar the electrodes, where is +-Volt/2.
# If set to 'ramp' is smoothly goes from -Volt/2 to +Volt/2 in a linear way.
#
print_log(" ")
print_log("    Initialization of V(i,j,k): "+how_to_initialize_V)
if how_to_initialize_V == 'ramp':
    initialize_V(Volt)
elif how_to_initialize_V == 'steps':
    set_electrodes(Volt)
else:
    set_electrodes(Volt)
#
# Preparing for parallel processing 5: converting the RHS of the equation from tensor to lists, to be used together with Veven, Vodd
# first of all, creating a list of indices corresponding to the boundary of the system, 
# namely i=0 or numofcellsalongx, j=0 or numofcellsalongy, k=0 or numofcellsalongz, on which the RHS is not evaluated.
#
print_log(" ")
evenboundaryindices=[]
for j in range(numofcellsalongy+1):
    for k in range(numofcellsalongz+1):
        num=0+j+k
        if num/2==num//2:
            evenboundaryindices.append(fromijktoevenoddindex[0][j][k])
for i in range(numofcellsalongx+1):
    for k in range(numofcellsalongz+1):
        num=i+0+k
        if num/2==num//2:
            evenboundaryindices.append(fromijktoevenoddindex[i][0][k])
for i in range(numofcellsalongx+1):
    for j in range(numofcellsalongy+1):
        num=i+j+0
        if num/2==num//2:
            evenboundaryindices.append(fromijktoevenoddindex[i][j][0])
for j in range(numofcellsalongy+1):
    for k in range(numofcellsalongz+1):
        num=numofcellsalongx+j+k
        if num/2==num//2:
            evenboundaryindices.append(fromijktoevenoddindex[numofcellsalongx][j][k])
for i in range(numofcellsalongx+1):
    for k in range(numofcellsalongz+1):
        num=i+numofcellsalongy+k
        if num/2==num//2:
            evenboundaryindices.append(fromijktoevenoddindex[i][numofcellsalongy][k])
for i in range(numofcellsalongx+1):
    for j in range(numofcellsalongy+1):
        num=i+j+numofcellsalongz
        if num/2==num//2:
            evenboundaryindices.append(fromijktoevenoddindex[i][j][numofcellsalongz])
oddboundaryindices=[]
for j in range(numofcellsalongy+1):
    for k in range(numofcellsalongz+1):
        num=0+j+k
        if num/2!=num//2:
            oddboundaryindices.append(fromijktoevenoddindex[0][j][k])
for i in range(numofcellsalongx+1):
    for k in range(numofcellsalongz+1):
        num=i+0+k
        if num/2!=num//2:
            oddboundaryindices.append(fromijktoevenoddindex[i][0][k])
for i in range(numofcellsalongx+1):
    for j in range(numofcellsalongy+1):
        num=i+j+0
        if num/2!=num//2:
            oddboundaryindices.append(fromijktoevenoddindex[i][j][0])
for j in range(numofcellsalongy+1):
    for k in range(numofcellsalongz+1):
        num=numofcellsalongx+j+k
        if num/2!=num//2:
            oddboundaryindices.append(fromijktoevenoddindex[numofcellsalongx][j][k])
for i in range(numofcellsalongx+1):
    for k in range(numofcellsalongz+1):
        num=i+numofcellsalongy+k
        if num/2!=num//2:
            oddboundaryindices.append(fromijktoevenoddindex[i][numofcellsalongy][k])
for i in range(numofcellsalongx+1):
    for j in range(numofcellsalongy+1):
        num=i+j+numofcellsalongz
        if num/2!=num//2:
            oddboundaryindices.append(fromijktoevenoddindex[i][j][numofcellsalongz])
#
# then defining the RHS in the list, rather than tensor, form.
#
def RHSeven(index):
    if evenboundaryindices.count(index)==0:
        return rhs(*fromevenoddindextoijk[0][index])
    else:
        return v(*fromevenoddindextoijk[0][index])
def RHSodd(index):
    if oddboundaryindices.count(index)==0:
        return rhs(*fromevenoddindextoijk[1][index])
    else:
        return v(*fromevenoddindextoijk[1][index])   
#
# Iterative solution of the equation, with output file updating
#
print_log('Calculation of R started...')
print_log(" ")
Rlist=[0]
#
# Estimate remaining time (measuring how long one single iteration takes)
#
start_iter_time = time.time()
for iteration in range(1):
        #sys.stdout.write('\r')
        #sys.stdout.write("    [%-50s] %d%%" % ('='*int(iteration*50/(iterations-1-0.001)), (5*20/(iterations-1-0.0001))*iteration)+"    Estimate: "+str(Rlist[-1]))
        #sys.stdout.flush()
            #
            # Calculate the even combination of the elements of V(i,j,k)
            #
        pool=mp.Pool(cpunumber)
        Veven=np.array(pool.map(RHSeven,range(dimVeven)))
        pool.close()
            #
        set_electrodes(Volt)
            #
            # Calculate the odd combination of the elements of V(i,j,k)
            #
        pool=mp.Pool(cpunumber)
        Vodd=np.array(pool.map(RHSodd,range(dimVodd)))
        pool.close()
            #
        set_electrodes(Volt)
            #
            # Compute R and append result to file
            #
        Rlist.append(calculate_R())
        with open(output_file_name_Resistance, "a+") as outputfile:
            outputfile.write("\n"+str(iteration)+","+str(Rlist[-1]))
        outputfile.close()
            #
stop_iter_time = time.time()
#
#
print_log('    This calculation will take about '+str(datetime.timedelta(seconds=round(iterations*(stop_iter_time-start_iter_time))))+' hh:mm:ss')
print_log('')
for iteration in range(iterations):
        sys.stdout.write('\r')
        sys.stdout.write("    [%-50s] %d%%" % ('='*int(iteration*50/(iterations-1-0.001)), (5*20/(iterations-1-0.0001))*iteration)+"    Estimate: "+str(Rlist[-1]))
        sys.stdout.flush()
            #
            # Calculate the even combination of the elements of V(i,j,k)
            #
        pool=mp.Pool(cpunumber)
        Veven=np.array(pool.map(RHSeven,range(dimVeven)))
        pool.close()
            #
        set_electrodes(Volt)
            #
            # Calculate the odd combination of the elements of V(i,j,k)
            #
        pool=mp.Pool(cpunumber)
        Vodd=np.array(pool.map(RHSodd,range(dimVodd)))
        pool.close()
            #
        set_electrodes(Volt)
            #
            # Compute R and append result to file
            #
        Rlist.append(calculate_R())
        with open(output_file_name_Resistance, "a+") as outputfile:
            outputfile.write("\n"+str(iteration)+","+str(Rlist[-1]))
        outputfile.close()
            #
print_log(" ")
print_log(" ")
print_log("Calculation completed.                  ")
#
# Output file so far: list of R for each iteration.
# Append to output file other outputs
#
#outputfile = open(outputfilename, 'a')
#outputfile.write('\n Cell phases and conductances, and nodes potential')
#outputfile.write('\n cell indices, cell phase, cellconductance along x, along y, along z, potential')
#outputfile.write('\n (i,j,k), phase(i,j,k), G_x(i,j,k), G_y(i,j,k), G_z(i,j,k), V(i,j,k)')
#for i in range(numofcellsalongx):
#    for j in range(numofcellsalongy):
#        for k in range(numofcellsalongz):
#            outputfile.write("\n("+str(i)+","+str(j)+","+str(k)+"),"+str(cellphase[i][j][k])+",("+str(gx[i,j,k])+","+str(gy[i,j,k])+","+str(gz[i,j,k])+"),"+str(v(i,j,k)))
#outputfile.close()
print_log(" ")
print_log("Estimated System Resistance: "+str(Rlist[-1])+" Ω.")
print_log(" ")
print_log("********************************************************************************************************")
#
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################)