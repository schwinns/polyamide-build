
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to replace charges')
args = parser.parse_args()

#######################################################################################
################################ STORE HEADER INFORMATION #############################
#######################################################################################

f = open(args.lmps, 'r')

# Box dimensions
for line in f:

    if line.startswith('Masses'):
        break

    elif 'xlo' in line.split(): # x box
        xlo = float(line.split()[0])
        xhi = float(line.split()[1])

    elif 'ylo' in line.split(): # y box
        ylo = float(line.split()[0])
        yhi = float(line.split()[1])

    elif 'zlo' in line.split(): # z box
        zlo = float(line.split()[0])
        zhi = float(line.split()[1])

# Atom masses
masses = {}
for line in f:

    if line.startswith('Atoms'):
        break

    elif len(line.split()) == 2:
        a_type = line.split()[0]
        mass = float(line.split()[1])

        masses[a_type] = mass


# Atoms section
total_mass = 0
zmin = 100
zmax = 0
for line in f:

    if line.startswith('Velocities'): # end of atoms section
        break

    elif line.startswith('Bonds'):
        break

    elif len(line.split()) > 0: # if not blank

        l = line.split()
        a_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        charge = l[3]
        x = float(l[4])
        y = float(l[5])
        z = float(l[6])

        if z < zmin:
            zmin = z
        if z > zmax:
            zmax = z

        total_mass += masses[a_type]


f.close()

#############################################################################################
####################################### QUICK CHECKS ########################################
#############################################################################################

total_mass = total_mass / 6.022 / 10**23 # [g/mol * mol/# = g]
s1 = (xhi-xlo) * 10**-8 # [Ang * 10^8 cm/Ang = cm]
s2 = (yhi-ylo) * 10**-8 # cm
s3 = (zhi-zlo) * 10**-8 # cm
vol = s1*s2*s3 # cm^3
density = total_mass / vol #/ 6.022 / 10

print('Density of the PA membrane using box dimensions: %.4f' %(density) )

s3 = (zmax-zmin) * 10**-8 # cm
vol = s1*s2*s3
density = total_mass / vol #/ 6.022 / 10

print('\nBox z dimension should be: %.4f %.4f zlo zhi' %(zmin, zmax) )
print('Density of the PA membrane using new dimensions: %.4f' %(density) )
