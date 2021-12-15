#!/usr/bin/env python

# Fix lmps file to have proper atom numbers after removing uncrosslinked oligomers

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='pdb file to rename atoms')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

# Define dictionary with possible bond, angle, dihedral types
# atom type 1, atom type 2 : bond type
bonds = {
    '1,1' : '1', # CA-CA
    '1,3' : '2', # CA-C
    '1,2' : '3', # CA-HA
    '3,4' : '4', # C-O
    '8,9' : '5', # LC-LN
    '1,5' : '6', # N-CA
    '5,6' : '7', # N-HN
    '1,7' : '8', # CA-NH
    '6,7' : '9', # NH-HN
    '1,8' : '10', # CA-LC
    '4,8' : '11', # LC-O
    '1,9' : '12', # LN-CA
    '6,9' : '13', # LN-HN
}


f = open(args.lmps,'r')
out = open(args.output, 'w')

# Header and Coeffs
for line in f:

    if line.startswith('Atoms'): # go to Atoms section
        out.write(line)
        break

    else: # everything before atoms
        out.write(line)

# Atoms section
atoms = {} # atom id : atom type
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        atom_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        charge = l[3]
        x = l[4]
        y = l[5]
        z = l[6]

        atoms[atom_id] = a_type

    else: # write blank lines
        out.write(line)

# Bonds section
for line in f:

    if line.startswith('Angles'): # go to Angles section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        b_id = l[0]
        old_type = l[1]
        a1_id = l[2]
        a2_id = l[3]
        a1_type = atoms[a1_id]
        a2_type = atoms[a2_id]

        if int(a1_type) <= int(a2_type):
            bond_pair = a1_type + ',' + a2_type
        else:
            bond_pair = a2_type + ',' + a1_type

        new_type = bonds[bond_pair]

        new_line = "%s %s %s %s\n" %(b_id, new_type, a1_id, a2_id)
        out.write(new_line)

    else: # write blank lines
        out.write(line)

# Angles section
for line in f:

    if line.startswith('Dihedrals'): # go to Dihedrals section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        old_type = l[1]
        a1_id = l[2]
        a2_id = l[3]
        a3_id = l[4]

        
    else: # write blank lines
        out.write(line)

# Dihedrals section
for line in f:

    if len(line.split()) > 0: # if line is not blank

        l = line.split()

        old_type = l[1]
        a1_id = l[2]
        a2_id = l[3]
        a3_id = l[4]
        a4_id = l[5]


    else: # write blank lines
        out.write(line)
