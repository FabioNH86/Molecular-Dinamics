##################################################################################################
# LUIS ADRIAN PADILLA SALAS
# THIS CODE GENERATES THE INITIAL CONFIGURATION FOR A SYSTEM WITH 2 SOLVENTS AND A GROUP OF 
# TRIBLOCK COPOLYMERS INITIALLY INTENDED FOR DPD SIMULATIONS.
# THIS CODE CREATES THE POLYMERS IN ONE SIDE OF THE SIMULATION BOX
# DEVELOPED BY MAY 19TH 2026
##################################################################################################

import random
import math
import numpy as np
from lxml import etree
import re
import sys
import numpy
numpy.set_printoptions(threshold=sys.maxsize)

# This code generates the initial configuration of a polymer immerse in a sea of solvent molecules
# such polymer has a sequence of A-B-C monomers where the number of monomers are NA - NB  respectively.
# thus the architecture is something like                A  A  A  A      A
#                                                        |  |  |  | .... |
#                                                        B  B  B  B .... B
#                                                        |  |  |  | .... |
#                                                        C  C  C  C .... C


# Parameters
Nc = 12000         # Number of chains in the system to be simulated
NA = 21            # Number of A monomers in an A branch
NB = 0            # Number of B monomers in a B branch
NC = 0            # Number of C monomers in a C branch
N = NA + NB + NC                    # Total number of monomers per polymer
Nmol = Nc * N                       # Total number of polymer beads 
n_bonds = Nmol - Nc                 # Total number of bonds on each brush polymer
rho = 3.0                           # Reduce density in DPD units
b = 0.75                            # Average bond length between connected polymer segments
massA = 1.0			    # Masses for the different DPD beads
massB = 1.0
massC = 1.0
massD = 1.0
massE = 1.0 

# In here we select if the system will create a subset of particles inside a droplet to help
# a fast thermalization of the system

droplet = 0
if droplet == 0:
   ############# Now the control is the polymer concentration
   Xp = 0.333                   # Polymer concentration
   Ntotal = int(Nmol*1.0/Xp)     # Total number of beads 
   Ntol = 1                   # Only one ceramic bead
   Naqua = Ntotal - Nmol - Ntol  # Total number of solvent beads  
   #############

# Simulation box parameters from selection of number of particles (elongated box)
ap = 4.0  #Aspect ratio of the simulation box Lx:Ly
ap2 = pow(ap, 2.0)
Volume = (1.0*Ntotal)/rho
LX3 = Volume*ap2 
Lx = pow(LX3, 1.0/3.0)
Ly = Lx/ap 
Lz = Lx/ap

print('Nmol', Nmol)
print('Ntol', Ntol)
print('Naqua', Naqua)
print('Ntotal', Ntotal)
print('Volume', Volume)
print('Lx', Lx)
print('Ly', Ly)
print('Lz', Lz)

# This data structure generates a unit vector
def random_three_vector():
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )

    x = x * b
    y = y * b
    z = z * b
    return (x,y,z)


# This data structure assigns the x, y and z coordinates for each bead
def positions():
    ab = NA + NB + NC   # Total number of A B molecules (in ONE Macro-molecule
    coor = np.zeros((ab, 3))

    # COORDINATES OF A, B & C BEADS
    ######################################################################################
    # Using the whole simulation box
    ######################################################################################
    if droplet == 0:
       for h in range(ab):
          x1, y1, z1 = random_three_vector()
          if h == 0:
             x = random.uniform(0.333, 1)
             x = x * (Lx/2)
             x = round(x, 4)
    
             y = random.uniform(-1, 1)
             y = y * (Ly/2)
             y = round(y, 4)
             
             z = random.uniform(-1, 1)
             z = z * (Lz/2)
             z = round(z, 4)
             
             coor[h, h] = x            
             coor[h, 1] = y
             coor[h, 2] = z
             
          else:
             x = x + x1
             x = x - Lx * int(round(x / Lx))
             x = round(x, 4)      

             y = y + y1
             y = y - Ly * int(round(y / Ly))
             y = round(y, 4)
             
             z = z + z1
             z = z - Lz * int(round(z / Lz))
             z = round(z, 4)
             
             coor[h, 0] = x
             coor[h, 1] = y
             coor[h, 2] = z
          
       return coor

# This data structure assigns the x, y and z coordinates for the solvent particles
def solvent_particles1():
    solv1 = np.zeros((Naqua, 3))

    #####################################################################################
    # Creating solvent over the whole simulation box
    #####################################################################################
    if droplet == 0:
       for p in range(Naqua):
          x = random.uniform(-1, 0.333)
          x = x * (Lx / 2)
          x = round(x, 4)
          
          y = random.uniform(-1, 1)
          y = y * (Ly / 2)
          y = round(y, 4)
          
          z = random.uniform(-1, 1)
          z = z * (Lz / 2)
          z = round(z, 4)
             
          solv1[p, 0] = x
          solv1[p, 1] = y
          solv1[p, 2] = z
    
       return solv1   
  
# This data structure assigns the x, y and z coordinates for the ceramic particles
def solvent_particles2():
    solv2 = np.zeros((Ntol, 3))

    #####################################################################################
    # Creating ceramic beads over the whole simulation box
    #####################################################################################
    if droplet == 0:
       for p in range(Ntol):
          x = random.uniform(-1, 1)
          x = x * (Lx / 2.0)
          x = round(x, 4)  

          y = random.uniform(-1, 1)
          y = y * (Ly / 2.0)
          y = round(y, 4)
           
          z = random.uniform(-1, 1)
          z = z * (Lz / 2.0)
          z = round(z, 4)
          
          solv2[p, 0] = x
          solv2[p, 1] = y
          solv2[p, 2] = z
          
       return solv2    

# This data structure assigns the molecule type to each position (1=A  2=B  3=C  D=solvent1 E=solvent2)
def type_assigning(x):
    rows = np.size(x, 0)
    print('arreglo', rows)

    y = np.chararray((rows, 1), unicode=True)
    ym = np.ones((rows, 1))
    y[:] = 'a'
    counter = 0
    beads = NA + NB + NC

    for i in range(Nc):

        for a in range(beads*counter, NA + beads*counter):
            y[a, 0] = 'A'
            ym[a, 0] = massA

        for b in range(NA + beads*counter, NA + NB + beads*counter):
            y[b, 0] = 'B'
            ym[b, 0] = massB

        for c in  range(NA + NB + beads*counter, NA + NB + NC + beads*counter):
            y[c, 0] = 'C'
            ym[c, 0] = massC

        counter = counter + 1
        #print(counter)

    count_aqua = 0
    for k in range(Nmol, Nmol+Naqua):
        y[k, 0] = 'D'
        ym[k, 0] = massD
        #count_aqua = count_aqua + 1
        #print(count_aqua, Naqua, Nmol)
 
    count_tol = 0
    for j in range(Nmol+Naqua, rows):
        y[j, 0] = 'E'
        ym[j, 0] = massE
        #count_tol = count_tol + 1
        #print(count_tol, Ntol)

    return y

# This data structure assigns the molecule mass to each position
def mass_assigning(x):
    rows = np.size(x, 0)

    ym = np.ones((rows, 1))
    counter = 0
    beads = NA + NB + NC

    for i in range(Nc):

        for a in range(beads*counter, NA + beads*counter):
            ym[a, 0] = massA

        for b in range(NA + beads*counter, NA + NB + beads*counter):
            ym[b, 0] = massB

        for c in  range(NA + NB + beads*counter, NA + NB + NC + beads*counter):
            ym[c, 0] = massC

        counter = counter + 1
        #print(counter)

    count_aqua = 0
    for k in range(Nmol, Nmol+Naqua):
        ym[k, 0] = massD
        #count_aqua = count_aqua + 1
        #print(count_aqua, Naqua, Nmol)
 
    count_tol = 0
    for j in range(Nmol+Naqua, rows):
        ym[j, 0] = massE
        #count_tol = count_tol + 1
        #print(count_tol, Ntol)

    return ym

def bonds_assigning():
    st = ''
    counter = 0
    ab = NA + NB + NC

    for l in range(Nc):
        check = 0

        for j in range(ab - 1):

            st += 'polymer' + '\t' + str(counter) + '\t' + str(counter + 1) + '\n'

            counter = counter + 1
            check = check + 1

        counter = counter + 1
        check = check + 1

    return st


def angle_assigning(y):
    string = ''
    counter = 0
    cc = 0

    for i in range(Nc):
        for j in range(N - 2):
            string += str(y[counter, 0]) + '-' + str(y[counter + 1, 0]) + '-' + str(y[counter + 2, 0]) + '\t'
            string += str(counter) + '\t' + str(counter + 1) + '\t' + str(counter + 2) + '\n'
            cc += 1
            counter = counter + 1
        counter += 2

    return string, cc

# if more than one Macro-molecule, use loop to concatenate positions calling
ww = positions()
for i in range(Nc - 1):
    wb = positions()
    ww = np.concatenate((ww, wb))

ss1 = solvent_particles1()
ss2 = solvent_particles2()

# concatenate Macro-molecules positions and solvent particles positions
waux = np.concatenate((ww, ss1))
w = np.concatenate((waux, ss2))

# assigns molecule type to the monomers
gh = type_assigning(w)

# get the number of columns and rows of the concatenated matrix
columns = np.size(w, 1)
rows = np.size(w, 0)

# This data structure creates an xml file using lmxl.etree module

# Create the root element
page = etree.Element('hoomd_xml', version='1.6')


# Add the page sub-element ' configuration'
configuration = etree.SubElement(page, 'configuration', natoms=str(Ntotal), dimensions='3', time_step='0')


# Add the configuration sub-elements
box = etree.SubElement(configuration, 'box', lz=str(Lz), ly=str(Ly), lx=str(Lx)).text = ''


position = etree.SubElement(configuration, 'position', num=str(Ntotal))
position.text = '\n' + re.sub("[][]", "", str(w)) + '\n'


velocity = etree.SubElement(configuration, 'velocity', num='0').text = ''


mass = etree.SubElement(configuration, 'mass', num=str(Ntotal))
m = mass_assigning(w)
#mass.text = '\n' + re.sub("[].[]", "", str(m)) + '\n'
#mass.text = '\n' + re.sub("[][.]", "", str(m)) + '\n'
mass.text = '\n' + re.sub("[][]", "", str(m)) + '\n'


diameter = etree.SubElement(configuration, 'diameter', num=str(Ntotal))
d = np.ones((rows, 1))
diameter.text = '\n' + re.sub('[][.]', "", str(d)) + '\n'


type = etree.SubElement(configuration, 'type', num=str(Ntotal))
type.text = '\n' + re.sub("[][']", "", str(gh)) + '\n'


body = etree.SubElement(configuration, 'body', num=str(Ntotal))
bo = np.ones((rows, 1))
for i in range(np.size(bo, 0)):
    bo[i] = -1
body.text = '\n' + re.sub("[][.]", "", str(bo)) + '\n'


bond = etree.SubElement(configuration, 'bond', num=str(n_bonds))
st = bonds_assigning()
bond.text = '\n' + st + '\n'

string, cc = angle_assigning(gh)
angle = etree.SubElement(configuration, 'angle', num=str(cc))
angle.text = '\n' + string + '\n'

dihedral = etree.SubElement(configuration, 'dihedral', num='0').text = ''

improper = etree.SubElement(configuration, 'improper', num='0').text = ''

charge = etree.SubElement(configuration, 'charge', num=str(Ntotal))
c = np.zeros((rows, 1))
charge.text = '\n' + re.sub("[][.]", "", str(c)) + '\n'

# Make a new document tree
doc = etree.ElementTree(page)

# Save to XML file
outFile = open('config.xml', 'wb')
doc.write(outFile, xml_declaration=True, pretty_print=True, encoding='utf-8')
