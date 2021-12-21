#!/usr/bin/env python

# Copy pdb coordinates to lmps file

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-p','--pdb',
                    help='pdb file with coordinates')
parser.add_argument('-l','--lmps',
                    help='lmps file with parameters')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

pdb = open(args.pdb, 'r')
lmp = open(args.lmps, 'r')
out = open(args.output, 'w')

# Write all the parameters and header from lammps file
for line in lmp:

    if line.startswith('Atoms'):
        out.write(line)
        out.write('\n')
        break

    else:
        out.write(line)

# Write the coordinates from pdb file in lammps format
n_waters = 1
waters = {}
n_atoms = 0
for line in pdb:

    if line.startswith('ATOM'): # write coordinates of water molecules

        l = line.split()
        a_id = l[1]
        a_type = l[2]
        mol_id = n_waters 
        x = l[6]
        y = l[7]
        z = l[8]

        n_atoms += 1

        # Assign types and charges for water molecules, create dictionary of water molecules
        if a_type == 'OW':
            waters[n_waters] = [n_atoms]
            a_type = '15'
            charge = '-0.830'
        elif a_type == 'HW':
            waters[n_waters].append(n_atoms)
            a_type = '16'
            charge = '0.415'
        else:
            print('Bad atom type: %s' %(a_type))
            exit()
        
        if len(waters[n_waters]) == 3:
            n_waters += 1

        new_line = " %s %s %s %s %s %s %s\n" %(n_atoms, mol_id, a_type, charge, x, y, z)
        out.write(new_line)
        
print('%d waters in the system' %(n_waters-1))
out.write('\nBonds\n\n')

# Write Bonds section
n_bonds = 0
for mol in waters:

    OW = waters[mol][0]

    for a in range(2):

        n_bonds += 1
        new_line = " %d %d %s %s\n" %(n_bonds, 19, OW, waters[mol][a+1])
        out.write(new_line)

out.write('\nAngles\n\n')

# Write Angles section
n_angles = 0
for mol in waters:

    n_angles += 1
    new_line = " %d %d %s %s %s\n" %(n_angles, 24, waters[mol][0], waters[mol][1], waters[mol][2])
    out.write(new_line)

print('\n%s atoms' %(n_atoms))
print('%d bonds' %(n_bonds))
print('%d angles' %(n_angles))