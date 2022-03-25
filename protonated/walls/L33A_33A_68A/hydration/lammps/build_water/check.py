#!/usr/bin/env python

# Add bottom reservoir of water

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lammps file')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

###############################################################################
############################ SAVE WATER TOPOLOGY ##############################
###############################################################################

lmp = open(args.lmps, 'r')

# Read in all information from lammps file
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
mols = {}
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
        mol_id = int(l[1])
        a_type = l[2]
        charge = l[3]
        x = float(l[4])
        y = float(l[5])
        z = float(l[6])

        if mol_id not in mols:
            mols[mol_id] = {
                # 'count'  : 1,
                'atoms'  : [a_id],
                'types'  : [a_type],
                'bonds'  : [],
                'angles' : []
            }
        else:
            # mols[mol_id]['count'] += 1
            mols[mol_id]['atoms'].append(a_id)
            mols[mol_id]['types'].append(a_type)

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

        mol_id = atoms[a1]['mol']
        mols[mol_id]['bonds'].append(b_id)

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

        mol_id = atoms[a1]['mol']
        mols[mol_id]['angles'].append(ang_id)

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

b = 0
a = 0
for m in sorted(mols):
    if not m == 1 and not len(mols[m]['bonds']) == 2:
        print('Molecule ',m, mols[m])
        b += 1
    if not m == 1 and not len(mols[m]['angles']) == 1:
        a += 1


print('\n%d waters do not have 2 bonds' %(b) )
print('%d waters do not have 1 angle' %(a) )