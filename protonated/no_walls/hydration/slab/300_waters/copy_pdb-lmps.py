#!/usr/bin/env python

# Copy pdb coordinates to lmps file

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-p','--pdb',
                    help='pdb file with coordinates')
parser.add_argument('-l','--lmps',
                    help='lmps file')
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

# Build dictionary with charges and atom ids --> charges[atom_id] = charge
charges = {}
for line in lmp:
    if line.startswith('Bonds'):
        vel = False
        break

    elif line.startswith('Velocities'):
        vel = True
        break

    elif len(line.split()) > 0:

        a_id = line.split()[0]
        charge = line.split()[3]

        charges[a_id] = charge

# Skip to bonding section
if vel:
    for line in lmp:
        if line.startswith('Bonds'):
            break

# Write the coordinates from pdb file in lammps format
n_waters = 1
waters = {}
maxes = [-100,-100,-100]
mins = [100,100,100]
n_atoms = 0
for line in pdb:

    if line.startswith('ATOM') and line.split()[3] == 'MOL': # write coordinates for the PA membrane and free OH's and H's

        l = line.split()
        a_id = l[1]
        a_type = l[2]
        mol_id = l[4]
        charge = charges[a_id]
        x = l[5]
        y = l[6]
        z = l[7]

        xyz = [float(x),float(y),float(z)]
        for d in range(len(xyz)):
            if xyz[d] > maxes[d]:
                maxes[d] = xyz[d]
            elif xyz[d] < mins[d]:
                mins[d] = xyz[d]

        n_atoms += 1

        new_line = " %d %s %s %s %s %s %s\n" %(n_atoms, mol_id, a_type, charge, x, y, z)
        out.write(new_line)
        n_mol = int(mol_id)

    elif line.startswith('ATOM') and len(line.split()) > 10: # write coordinates of water molecules

        l = line.split()
        a_id = l[1]
        a_type = l[2]
        mol_id = n_mol + n_waters # using last mol_id from PA section
        x = l[6]
        y = l[7]
        z = l[8]

        xyz = [float(x),float(y),float(z)]
        for d in range(len(xyz)):
            if xyz[d] > maxes[d]:
                maxes[d] = xyz[d]
            elif xyz[d] < mins[d]:
                mins[d] = xyz[d]

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

    elif line.startswith('ATOM'): # write coordinates of water molecules where molecule number touches A/B

        l = line.split()
        a_id = l[1]
        a_type = l[2]
        mol_id = n_mol + n_waters # using last mol_id from PA section
        x = l[5]
        y = l[6]
        z = l[7]

        xyz = [float(x),float(y),float(z)]
        for d in range(len(xyz)):
            if xyz[d] > maxes[d]:
                maxes[d] = xyz[d]
            elif xyz[d] < mins[d]:
                mins[d] = xyz[d]

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
        
print('%d waters in the system\n' %(n_waters-1))
print('Box dimensions should be:')
print('\t%f %f xlo xhi' %(mins[0], maxes[0]))
print('\t%f %f ylo yhi' %(mins[1], maxes[1]))
print('\t%f %f zlo zhi' %(mins[2], maxes[2]))

# Write Bonds section
out.write('\nBonds\n\n')
for line in lmp:

    if line.startswith('Angles'):
        break

    elif len(line.split()) > 0:
        n_bonds = int(line.split()[0])
        out.write(line)


# Add new bonds for waters
for mol in waters:

    OW = waters[mol][0]

    for a in range(2):

        n_bonds += 1
        new_line = " %d %d %s %s\n" %(n_bonds, 19, OW, waters[mol][a+1])
        out.write(new_line)

out.write('\nAngles\n\n')

# Write Angles section
for line in lmp:

    if line.startswith('Dihedrals'):
        break

    elif len(line.split()) > 0:
        n_angles = int(line.split()[0])
        out.write(line)


# Add new angles for waters
for mol in waters:

    n_angles += 1
    new_line = " %d %d %s %s %s\n" %(n_angles, 24, waters[mol][0], waters[mol][1], waters[mol][2])
    out.write(new_line)

out.write('\nDihedrals\n')

# Write Dihedrals section
for line in lmp:
    out.write(line)

print('\n%s atoms' %(n_atoms))
print('%d bonds' %(n_bonds))
print('%d angles' %(n_angles))