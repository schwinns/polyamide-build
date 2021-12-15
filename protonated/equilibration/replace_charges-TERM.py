#!/usr/bin/env python

# Fix lmps file to have proper atom charges after polymerization
#   Assuming there are args.MPD MPD molecules, followed by args.TMC_L TMC-L molecules,
#   followed by args.TMC_C TMC_C molecules

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to replace charges')
parser.add_argument('-m','--MPD',default=247,
                    help='number of MPD molecules in system')
parser.add_argument('-t','--TMC_L',default=101,
                    help='number of MPD molecules in system')
parser.add_argument('-c','--TMC_C',default=97,
                    help='number of MPD molecules in system')                    
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

# Define atom charges as those from GAFF reparameterization
#       charges[monomer][index] = charge
charges = {
    'MPD_L' : { # comments are names in 1MPD-2TMC to show fully crosslinked charges
        1  : -0.1875, # C7
        2  : -0.079, # C9
        3  : -0.1875, # C10
        4  : 0.0896, # C8
        5  : -0.206, # C12
        6  : 0.0896, # C11
        7  : 0.1545, # H3
        8  : 0.14, # H4
        9  : 0.1545, # H5
        10 : -0.4651, # N1
        11 : -0.4651, # N
        12 : 0.17, # H6 
        13 : 0.324, # H7
        14 : 0.324 # H10
    },
    'TMC_L' : { # comments are names in 2MPD-1TMC to show fully crosslinked charges
        1  : -0.031, # C2
        2  : -0.1486, # C
        3  : -0.031, # C1
        4  : -0.1686, # C3
        5  : 0.003, # C4
        6  : -0.1686, # C5
        7  : 0.167, # H1
        8  : 0.167, # H
        9  : 0.184, # H2
        10 : 0.6802, # C6
        11 : 0.6802, # C14
        12 : 0.6507, # C13
        13 : -0.543, # O2
        14 : -0.5706, # O
        15 : -0.5706, # O1
        16 : -0.6071, # O3
        17 : 0.453 # H17
    },
    'TMC_C' : { # comments are names in 3MPD-1TMC to show fully crosslinked charges
        1  : -0.038333, # C4
        2  : -0.159267, # C5
        3  : -0.038333, # C1
        4  : -0.159267, # C3
        5  : -0.038333, # C2
        6  : -0.159267, # C
        7  : 0.165333, # H2
        8  : 0.165333, # H
        9  : 0.165333, # H1
        10 : 0.679367, # C6
        11 : 0.679367, # C13
        12 : 0.679367, # C14
        13 : -0.5691, # O1
        14 : -0.5691, # O
        15 : -0.5691 # O2
    },
    'term' : {
        '1m' : 0.1686, # CA-MPD
        1    : -0.1581, # CA-TMC
        4    : -0.5425, # O
        6    : 0.3948, # HN
        7    : -0.6036, # OH
        8    : -0.8232, # NH
        9    : 0.453, # HO
        12   : 0.6517  # CT
    } 
}

atoms_per_MPD = len(charges['MPD'])
atoms_per_TMC_L = len(charges['TMC_L'])
atoms_per_TMC_C = len(charges['TMC_C'])

# Keep track of the total charge in the system and how many of each atom type
total_charge = 0
n_types = {
    '1'  : 0, # CA
    '2'  : 0, # HA
    '3'  : 0, # C
    '4'  : 0, # O
    '5'  : 0, # N
    '6'  : 0, # HN
    '7'  : 0, # OH
    '8'  : 0, # NH
    '9'  : 0, # HO
    '10' : 0, # LC
    '11' : 0, # LN
    '12' : 0, # CT
}

# Define the atom types associated with terminated groups
terminated = ['7','8','12']

# Define bond types to ignore when saving atoms
ignore_bonds = ['8','16'] # CA-NH, CA-CT

# Find all the terminated groups and explicitly set their charges
term_atoms = []
f = open(args.lmps, 'r')
for line in f:
    if line.startswith('Atoms'):
        break

for line in f: # find all the terminated atom types and save their atom ids
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

        if a_type in terminated:
            term_atoms.append(a_id)
    
if vel:
    for line in f:
        if line.startswith('Bonds'):
            break

for line in f: # find the atom ids bonded to terminated atoms
    
    if line.startswith('Angles'):
        break

    elif len(line.split()) > 0:

        l = line.split()
        b_id = l[0]
        b_type = l[1]
        a1 = l[2]
        a2 = l[3]

        if a1 in term_atoms and a2 not in term_atoms and b_type not in ignore_bonds:
            term_atoms.append(a2)
        elif a2 in term_atoms and a1 not in term_atoms and b_type not in ignore_bonds:
            term_atoms.append(a1)
        
f.close()

f = open(args.lmps,'r')
out = open(args.output, 'w')

# Write Header and Coeffs
for line in f:

    if line.startswith('Atoms'): # go to Atoms section
        out.write(line)
        break

    else: # everything before atoms
        out.write(line)

# Atoms section
i = 1
stop = False
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        a_id = int(l[0])
        mol_id = l[1]
        a_type = l[2]

        if str(a_id) in term_atoms: # if we are a termination atom
            monomer = 'term'
            idx = int(a_type)
            if idx == 1 and i <= int(args.MPD) * atoms_per_MPD: # aromatic carbon bonded to terminated group is exception (MPD is different than TMC)
                idx = '1m'
            i += 1
        elif i <= int(args.MPD) * atoms_per_MPD: # if we are in the first args.MPD*atoms_per_MPD 
            monomer = 'MPD'
            idx = i % atoms_per_MPD
            if idx == 0:
                idx = int(atoms_per_MPD)
            i += 1
        elif i <= int(args.TMC_L) * atoms_per_TMC_L + int(args.MPD) * atoms_per_MPD: # if we are in the TMC_L monomers
            monomer = 'TMC_L'
            idx = (i - int(args.MPD) * atoms_per_MPD) % atoms_per_TMC_L
            if idx == 0:
                idx = int(atoms_per_TMC_L)
            i += 1
        else:
            monomer = 'TMC_C'
            idx = (i - int(args.MPD) * atoms_per_MPD - int(args.TMC_L) * atoms_per_TMC_L) % atoms_per_TMC_C
            if idx == 0:
                idx = int(atoms_per_TMC_C)
            i += 1


        charge = charges[monomer][idx]
        total_charge += charge
        n_types[a_type] += 1
        
        print('ID: %d, index: %s' %(a_id,str(idx)))
        print(charge)
        # print('Total charge: %.4f' %(total_charge))
        if stop:
            exit()

        x = l[4]
        y = l[5]
        z = l[6]

        new_line = " %s %s %s %.4f %s %s %s\n" %(a_id,mol_id,a_type,charge,x,y,z)
        out.write(new_line)

    else: # write blank lines
        out.write(line)

print('Total charge in the system is %.4f' %(total_charge))
print('The number of each atom type:')
print(n_types)

# Write all other sections unchanged
for line in f:
    out.write(line)
