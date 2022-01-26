#!/usr/bin/env python

# Combine lmps files with water and PA membrane

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-pa','--lmps_PA',
                    help='lammps file with PA membrane')
parser.add_argument('-w','--lmps_water',
                    help='lmps file with water box')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

###############################################################################
############################ SAVE WATER TOPOLOGY ##############################
###############################################################################

lmp = open(args.lmps_water, 'r')

# Read in all information from water lammps file to add to end of PA membrane file
masses = {}
pair_coeffs = {}
for line in lmp:

    if line.startswith('Atoms'):
        break

    elif len(line.split()) == 2:
        # Number of atoms, bonds, etc
        if 'atoms' in line.split():
            n_atoms = int(line.split()[0])
        elif 'bonds' in line.split():
            n_bonds = int(line.split()[0])
        elif 'angles' in line.split():
            n_angles = int(line.split()[0])
        elif 'dihedrals' in line.split():
            n_dihedrals = int(line.split()[0])
        
        # Masses
        else:
            masses[line.split()[0]] = line.split()[1]

    elif len(line.split()) == 3:
        # Types of atoms, bonds, etc
        if 'atom' in line.split():
            a_types = int(line.split()[0])
        elif 'bond' in line.split():
            b_types = int(line.split()[0])
        elif 'angle' in line.split():
            ang_types = int(line.split()[0])
        elif 'dihedral' in line.split():
            dih_types = int(line.split()[0])

        # Pair Coeffs
        else:
            pair_coeffs[line.split()[0]] = line.split()[1:]

    elif len(line.split()) == 4 and '#' not in line.split(): # box dimensions
        if 'xlo' in line.split():
            xlo = float(line.split()[0])
            xhi = float(line.split()[1])
        elif 'ylo' in line.split():
            ylo = float(line.split()[0])
            yhi = float(line.split()[1])
        elif 'zlo' in line.split():
            zlo = float(line.split()[0])
            zhi = float(line.split()[1])

# Build dictionary with atomistic topology -- atoms section
atoms = {}
for line in lmp:
    if line.startswith('Bonds'):
        vel = False
        break

    elif line.startswith('Velocities'):
        vel = True
        break

    elif len(line.split()) > 0:
        
        l = line.split()
        a_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        charge = l[3]
        x = float(l[4])
        y = float(l[5])
        z = float(l[6])

        atoms[a_id] = {
            'mol'    : mol_id,
            'type'   : a_type,
            'charge' : charge,
            'coords' : np.array([x,y,z]),
            'bonded' : []
        }

# Skip to bonding section
if vel:
    for line in lmp:
        if line.startswith('Bonds'):
            break

# Build dictionary with atomistic topology -- bonds section
bonds = {}
for line in lmp:

    if line.startswith('Angles'):
        break

    elif len(line.split()) > 0:

        l = line.split()
        b_id = l[0]
        b_type = l[1]
        a1 = l[2]
        a2 = l[3]

        atoms[a1]['bonded'].append(a2)
        atoms[a2]['bonded'].append(a1)

        bonds[b_id] = {
            'type'  : b_type,
            'atoms' : [a1,a2]
        }

# Build dictionary with atomistic topology -- angles section
angles = {}
for line in lmp:

    if line.startswith('Dihedrals'):
        break

    elif len(line.split()) > 0:

        l = line.split()
        ang_id = l[0]
        ang_type = l[1]
        a1 = l[2]
        a2 = l[3]
        a3 = l[4]

        angles[ang_id] = {
            'type'  : ang_type,
            'atoms' : [a1,a2,a3]
        }

# Build dictionary with atomistic topology -- dihedrals section
dihedrals = {}
for line in lmp:

    if line.startswith('Impropers'):
        break

    elif len(line.split()) > 0:

        l = line.split()
        dih_id = l[0]
        dih_type = l[1]
        a1 = l[2]
        a2 = l[3]
        a3 = l[4]
        a4 = l[5]

        dihedrals[dih_id] = {
            'type'  : dih_type,
            'atoms' : [a1,a2,a3,a4]
        }

if len(atoms) != n_atoms:
    print('Numbers of atoms do not match...')
    exit()
if len(bonds) != n_bonds:
    print('Numbers of bonds do not match...')
    exit()
if len(angles) != n_angles:
    print('Numbers of angles do not match...')
    exit()
if len(dihedrals) != n_dihedrals:
    print('Numbers of dihedrals do not match...')
    exit()

###############################################################################
######################### ADD WATERS TO PA MEMBRANE ###########################
###############################################################################

out = open(args.output, 'w')
lmp = open(args.lmps_PA, 'r')

# Write header
m = 0
for line in lmp:

    if line.startswith('Atoms'):
        out.write(line)
        break

    elif len(line.split()) == 2:
        # Number of atoms, bonds, etc
        if 'atoms' in line.split():
            n = n_atoms + int(line.split()[0])
            out.write('%d atoms\n' %(n))
        elif 'bonds' in line.split():
            n = n_bonds + int(line.split()[0])
            out.write('%d bonds\n' %(n))
        elif 'angles' in line.split():
            n = n_angles + int(line.split()[0])
            out.write('%d angles\n' %(n))
        elif 'dihedrals' in line.split():
            n = n_dihedrals + int(line.split()[0])
            out.write('%d dihedrals\n' %(n))
        
        # Masses
        else:
            m += 1
            mass = masses[str(m)]
            out.write('%d %s\n' %(m, masses[str(m)]) )

    elif len(line.split()) == 3:
        # Types of atoms, bonds, etc
        if 'atom' in line.split():
            out.write('%s atom types\n' %(a_types) )
        elif 'bond' in line.split():
            out.write('%s bond types\n' %(b_types) )
        elif 'angle' in line.split():
            out.write('%s angle types\n' %(ang_types) )
        elif 'dihedral' in line.split():
            out.write('%s dihedral types\n' %(dih_types) )

        # Pair Coeffs
        else:
            pair_coeffs[line.split()[0]] = line.split()[1:]

    elif len(line.split()) == 4 and '#' not in line.split(): # box dimensions
        if 'xlo' in line.split():
            lo = float(line.split()[0])
            hi = float(line.split()[1])
            cross_x = hi - lo
            out.write(line)
        elif 'ylo' in line.split():
            lo = float(line.split()[0])
            hi = float(line.split()[1])
            cross_y = hi - lo
            out.write(line)
        elif 'zlo' in line.split():
            zlo = float(line.split()[0])
            zhi = float(line.split()[1])

    else:
        out.write(line)

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
        new_line = "%d %d %s %s\n" %(n_bonds, 19, OW, waters[mol][a+1])
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
    new_line = "%d %d %s %s %s\n" %(n_angles, 24, waters[mol][0], waters[mol][1], waters[mol][2])
    out.write(new_line)

out.write('\nDihedrals\n')

# Write Dihedrals section
for line in lmp:
    out.write(line)

print('\n%s atoms' %(n_atoms))
print('%d bonds' %(n_bonds))
print('%d angles' %(n_angles))