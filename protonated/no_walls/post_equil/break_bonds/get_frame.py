# Fix lmps file to have parameters and bonding information

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-hd','--header',
                    help='lmps file with header and box dimensions')
parser.add_argument('-p','--params',
                    help='lmps file with parameters')
parser.add_argument('-d','--data',
                    help='lmps file with coordinate data')
parser.add_argument('-b','--bonds',
                    help='lmps file with the bonding information')
parser.add_argument('-t','--timestep',default=-1,type=int,
                    help='timestep to pull from the lammps trajectory, default is last frame')
parser.add_argument('-bi','--box_image',action='store_true',
                    help='output writes the box image of each atom')
parser.add_argument('-v','--velocity',action='store_true',
                    help='output writes the velocities of each atom')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

#######################################################################################
############################### WRITE THE HEADER INFORMATION ##########################
#######################################################################################

head = open(args.header,'r')
out = open(args.output, 'w')

# Header
for line in head:

    if line.startswith('Masses'): # go to Parameters section
        out.write(line)
        break

    else: # everything before masses
        out.write(line)

out.write('\n')
head.close()
params = open(args.params, 'r')

#######################################################################################
############################ WRITE PARAMETERS AND SAVE MASSES #########################
#######################################################################################

# Parameters section
for line in params: # skip until Masses
    if line.startswith('Masses'):
        break

masses = {} # save the masses for density calculation
for line in params:

    if line.startswith('Atoms'): # go to Atoms section
        out.write(line)
        break

    elif len((line.split('#')[0]).split()) == 2 and 'Coeffs' not in line: # if it is a mass
        out.write(line)
        l = line.split()
        masses[l[0]] = float(l[1])

    else: # everything before masses
        out.write(line)

out.write('\n')
params.close()
data = open(args.data, 'r')

#######################################################################################
########################## FIND THE DATA FOR THE CHOSEN TIMESTEP ######################
#######################################################################################

# Atoms section
check = False
ts = 0
for line in data: # find the chosen timestep

    if line.startswith('ITEM: TIMESTEP'):
        check = True
    
    elif check:
        ts = int(line.split()[0])
        check = False

    if ts == args.timestep:
        break
    elif args.timestep == -1:
        last_ts = ts

if args.timestep == -1: # find the last saved timestep

    print('Final timestep is %d' %(last_ts))

    data.close()
    data = open(args.data, 'r')

    check = False
    ts = 0
    for line in data: # skip everything until the last timestep
        if line.startswith('ITEM: TIMESTEP'):
            check = True

        elif check:
            ts = int(line.split()[0])
            check = False

        if ts == last_ts:
            break

#######################################################################################
################################## SAVE BOX DIMENSIONS ################################
#######################################################################################

dims = []
for line in data:
    if line.startswith('ITEM: BOX'):
        break

for line in data: # skip info before atom coordinates

    if line.startswith('ITEM: ATOMS'): # go to Bonds section
        break

    else:
        l = line.split()
        lo = float(l[0])
        hi = float(l[1])

        dims.append([lo,hi])

print('Box dimensions (Angstroms):')
bl = []
for dim in dims:
    print('\t%f %f' %(dim[0], dim[1]))
    bl.append(dim[1] - dim[0])

vol = (bl[0] * 10**-8)**3
print('Volume (cm^3):')
print('\t%e' %(vol))

#######################################################################################
##################################### WRITE ATOM DATA #################################
#######################################################################################

out.write('\nAtoms\n\n')
mass = 0
for line in data: # write all the atom coordinates
    
    l = line.split()
    a_id = l[0]
    mol_id = l[1]
    a_type = l[2]
    charge = l[3]
    xyz = l[4] + ' ' + l[5] + ' ' + l[6]
    vel = l[7] + ' ' + l[8] + ' ' + l[9]
    box_im = l[10] + ' ' + l[11] + ' ' + l[12]

    mass += masses[a_type]

    if args.velocity and args.box_image:
        new_line = line
    elif args.velocity:
        new_line = '%s %s %s %s %s %s\n' %(a_id, mol_id, a_type, charge, xyz, vel)
    elif args.box_image:
        new_line = '%s %s %s %s %s %s\n' %(a_id, mol_id, a_type, charge, xyz, box_im)
    else:
        new_line = '%s %s %s %s %s\n' %(a_id, mol_id, a_type, charge, xyz)

    out.write(new_line)

mass = mass / 6.022 / 10**23
print('Mass (g):')
print('\t%e' %(mass))
density = mass / vol
print('Density (g/cm^3):')
print('\t%e\n' %(density))

#######################################################################################
################################ WRITE BONDING DATA ###################################
#######################################################################################

out.write('\n')
data.close()
bonds = open(args.bonds, 'r')

# Bonds, angle, dihedral information
for line in bonds: # skip everything before the bonding info

    if line.startswith('Bonds'):
        out.write(line)
        break

for line in bonds: # write all the bond information
    out.write(line)

