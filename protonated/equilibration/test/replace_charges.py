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

#############################################################################################
####################################### USER INPUTS #########################################
#############################################################################################

# Define atom charges as those from GAFF reparameterization
#       charges[monomer][index] = charge
charges = {
    'MPD_L' : { # comments are names in 1MPD-2TMC
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
    'TMC_L' : { # comments are names in 2MPD-1TMC
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
    'TMC_C' : { # comments are names in 3MPD-1TMC
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
    'MPD_E' : { # comments are names in 1MPD-1TMC
        1  : -0.223, # C7
        2  : -0.067, # C9
        3  : -0.215, # C10
        4  : 0.0996, # C8
        5  : -0.216, # C12
        6  : 0.1686, # C11
        7  : 0.13, # H3
        8  : 0.131, # H4
        9  : 0.134, # H5
        10 : -0.8232, # N1
        11 : -0.4651, # N
        12 : 0.17, # H6
        13 : 0.3215, # H7
        14 : 0.3948 # H8
    },
    'TMC_E' : { # comments are names in 1MPD-1TMC
        1  : 0.008, # C1
        2  : -0.1581, # C5
        3  : -0.006, # C4
        4  : -0.1581, # C
        5  : -0.006, # C2
        6  : -0.1766, # C3
        7  : 0.184, # H
        8  : 0.178, # H2
        9  : 0.178, # H1
        10 : 0.6517, # C13
        11 : 0.6807, # C6
        12 : 0.6517, # C14
        13 : -0.5425, # O1
        14 : -0.5425, # O3
        15 : -0.5761, # O
        16 : -0.6036, # O2
        17 : 0.453 # H11
    },
    'term' : {
        6    : 0.3948, # HN
        7    : -0.6036, # OH
        9    : 0.453 # HO
    } 
}

atoms_per_MPD_L = len(charges['MPD_L'])
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
terminated = ['8','12']

#############################################################################################
############################### CLASSIFICATION OF MONOMERS ##################################
#############################################################################################

# Find all the terminated groups and classify each monomer
f = open(args.lmps, 'r')
for line in f:
    if line.startswith('Atoms'):
        break

monomers = {}
mono = 1
i = 1
for line in f:

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

        if i <= int(args.MPD) * atoms_per_MPD_L: # if we are in the first args.MPD*atoms_per_MPD_L 
            
            if i % atoms_per_MPD_L == 1: # if first atom in monomer make new entry in monomers dict
                monomers[mono] = {
                    'atoms'  : [], # save atoms
                    'type'   : 'MPD_L',  # save type of monomer
                    'n_term' : 0 # save number of terminated N's
                }

            if a_type in terminated: # change monomer type if a terminated atom is included
                monomers[mono]['type'] = 'MPD_E'
                monomers[mono]['n_term'] += 1

            monomers[mono]['atoms'].append(a_id) # add the atom id to the atoms list

            if i % atoms_per_MPD_L == 0: # if last atom in monomer, increment mono
                mono += 1

            i += 1

        elif i <= int(args.TMC_L) * atoms_per_TMC_L + int(args.MPD) * atoms_per_MPD_L: # if we are in the TMC_L monomers
            
            if (i - int(args.MPD) * atoms_per_MPD_L) % atoms_per_TMC_L == 1: # if first atom in monomer make new entry in monomers dict
                monomers[mono] = {
                    'atoms' : [],
                    'type'  : 'TMC_L',
                    'n_term' : 0
                }
                
            if a_type in terminated:
                monomers[mono]['n_term'] += 1
            
            if monomers[mono]['n_term'] == 2: # if two of reactive C's are terminated, change monomer type
                monomers[mono]['type'] = 'TMC_E'

            monomers[mono]['atoms'].append(a_id) # add the atom id to the atoms list

            if (i - int(args.MPD) * atoms_per_MPD_L) % atoms_per_TMC_L == 0:
                mono += 1

            i += 1

        elif i <= int(args.MPD) * atoms_per_MPD_L + int(args.TMC_L) * atoms_per_TMC_L + int(args.TMC_C) * atoms_per_TMC_C:

            if (i - int(args.MPD) * atoms_per_MPD_L - int(args.TMC_L) * atoms_per_TMC_L) % atoms_per_TMC_C == 1: # if first atom in monomer make new entry in monomers dict
                monomers[mono] = {
                    'atoms' : [],
                    'type'  : 'TMC_C',
                    'n_term' : 0
                }
                
            if a_type in terminated:
                monomers[mono]['n_term'] += 1
            
            if monomers[mono]['n_term'] == 1: # if one of reactive C's are terminated, change monomer type to TMC_L
                monomers[mono]['type'] = 'TMC_L'
            elif monomers[mono]['n_term'] == 2: # if two of reactive C's are terminated, change monomer type to TMC_E
                monomers[mono]['type'] = 'TMC_C'

            monomers[mono]['atoms'].append(a_id) # add the atom id to the atoms list

            if (i - int(args.MPD) * atoms_per_MPD_L - int(args.TMC_L) * atoms_per_TMC_L) % atoms_per_TMC_C == 0:
                mono += 1

            i += 1

        else: # if it is a termination atom added to the end of the Atoms section

            if (i - int(args.MPD) * atoms_per_MPD_L - int(args.TMC_L) * atoms_per_TMC_L) % atoms_per_TMC_C == 1: # if first atom in monomer make new entry in monomers dict
                monomers[mono] = {
                    'atoms' : [],
                    'type'  : 'term',
                    'n_term' : 0
                }

            monomers[mono]['atoms'].append(a_id)

            i += 1

print(monomers)
f.close()

#############################################################################################
################################### REWRITE WITH CHARGES ####################################
#############################################################################################

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
mono = 1
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        vel = False
        break

    elif line.startswith('Velocities'):
        vel = True
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        a_id = l[0]
        mol_id = l[1]
        a_type = l[2]

        if a_id in monomers[mono]['atoms']: # if atom id is in list of atoms in monomer
            monomer = monomers[mono]['type']
            idx = monomers[mono]['atoms'].index(a_id) + 1

        else: # otherwise check the next monomer
            mono += 1

        if a_id in monomers[mono]['atoms']:
            monomer = monomers[mono]['type']
            idx = monomers[mono]['atoms'].index(a_id) + 1

        if monomers[mono]['type'] == 'term': # if the atom is an added termination atom
            monomer = monomers[mono]['type']
            idx = int(a_type)

        charge = charges[monomer][idx]
        total_charge += charge
        n_types[a_type] += 1

        # print(a_id)
        # print(monomer)
        # print(idx)
        # print(charge)
        # print(total_charge)
        # print()
        
        x = l[4]
        y = l[5]
        z = l[6]

        new_line = " %s %s %s %.4f %s %s %s\n" %(a_id,mol_id,a_type,charge,x,y,z)
        out.write(new_line)

    else: # write blank lines
        out.write(line)

# Skip velocities if present
if vel:
    for line in f:
        if line.startswith('Bonds'):
            break

# Write all other sections unchanged
for line in f:
    out.write(line)

#############################################################################################
####################################### QUICK CHECKS ########################################
#############################################################################################

print('Total charge in the system is %.4f' %(total_charge))
print('The number of each atom type:')
print(n_types)