#!/usr/bin/env python

# Fix lmps file to have proper atom numbers for packing
# Must manually change atom types in lmps Atoms section

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='pdb file to rename atoms')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

# Define dictionaries with possible bond, angle, dihedral types
# atom type 1, atom type 2 : bond type (atom types must be in increasing order)
bonds = {
    '1,1'   : '1', # CA-CA
    '1,3'   : '2', # CA-C
    '1,2'   : '3', # CA-HA
    '3,4'   : '4', # C-O
    '10,11' : '5', # LC-LN
    '1,5'   : '6', # N-CA
    '5,6'   : '7', # N-HN
    '1,8'   : '8', # CA-NH
    '6,8'   : '9', # NH-HN
    '1,10'  : '10', # CA-LC
    '4,10'  : '11', # LC-O
    '1,11'  : '12', # LN-CA
    '6,11'  : '13', # LN-HN
    '7,12'  : '14', # CT-OH
    '7,9'   : '15', # OH-HO
    '1,12'  : '16', # CA-CT
    '4,12'  : '17'  # CT-O
}

# atom type 1, atom type 2, atom type 3 : angle type
angles = {
    '1,1,1' : '1', # CA-CA-CA
    '1,1,2' : '2', # CA-CA-HA
    '1,3,4' : '3', # CA-C-O
    '1,1,3' : '4', # CA-CA-C
    '1,10,11' : '5', # CA-LC-LN (these two are indistinguishable in the code)
    '1,10,11' : '6', # LC-LN-CA (but it does not matter because they will not appear in monomers)
    '6,10,11' : '7', # LC-LN-HN
    '4,10,11' : '8', # O-LC-LN
    '1,1,5' : '9', # N-CA-CA
    '1,5,6' : '10', # CA-N-HN
    '1,1,8' : '11', # CA-CA-NH
    '1,6,8' : '12', # CA-NH-HN
    '3,4,4' : '13', # O-C-O
    '6,6,8' : '14', # HN-NH-HN
    '1,4,10' : '15', # CA-LC-O
    '1,1,10' : '16', # CA-CA-LC
    '1,1,11' : '17', # LN-CA-CA
    '1,6,11' : '18', # CA-LN-HN
    '1,7,12' : '19', # CA-CT-OH
    '7,9,12' : '20', # CT-OH-HO
    '4,7,12' : '21', # O-CT-OH
    '1,4,12' : '22', # CA-CT-O
    '1,1,12' : '23'  # CA-CA-Ct
}

# atom types : dih type
dihedrals = {
    '1,1,1,1' : '1',  # CA-CA-CA-CA
    '1,1,1,2' : 'check1',  # CA-CA-CA-HA (improper)
    '1,1,3,4' : '3',  # CA-CA-C-O
    '1,1,10,11' : '4',  # CA-LC-LN-CA (same as 6,30, but will not appear)
    '1,4,10,11' : '5',  # CA-LC-LN-O (same as 32, but will not appear)
    '1,1,10,11' : '6',  # LC-LN-CA-CA (same as 4,30)
    '4,6,10,11' : '7',  # O-LC-LN-HN (same as 33, but will not appear)
    '1,1,6,8' : '8',  # CA-CA-NH-HN
    '1,1,1,3' : 'check2',  # CA-CA-CA-C 
    '1,1,1,10' : '10', # CA-CA-CA-LC (same as 23, but will not appear)
    # '1,1,1,2' : '11', # CA-CA-CA-HA (same as check1)
    '1,1,2,3' : '12', # C-CA-CA-HA
    '1,1,2,10' : '13', # LC-CA-CA-HA
    '1,1,1,5' : 'check3', # N-CA-CA-CA
    '1,1,1,11' : '15', # LN-CA-CA-CA (same as 26, but will not appear)
    '1,1,2,5' : '16', # N-CA-CA-HA
    '1,1,2,11' : '17', # LN-CA-CA-HA
    '1,1,1,8' : 'check4', # CA-CA-CA-NH
    '1,1,2,2' : '19', # HA-CA-CA-HA
    '1,1,2,8' : '20', # HA-CA-CA-NH
    '1,4,7,12' : '21', # CA-O-OH-CT (improper)
    # '1,1,1,3' : '22',  # CA-CA-CA-C (improper) (same as check2)
    '1,1,1,10' : '23',  # CA-CA-CA-LC (improper) (same as 10)
    '1,6,10,11' : 'check5',  # LC-CA-LN-HN (improper)
    # '1,1,1,5' : '25',  # CA-CA-CA-N (improper) (same as check3)
    '1,1,1,11' : '26',  # CA-CA-CA-LN (improper) (same as 15)
    # '1,1,1,8' : '27',  # CA-CA-CA-NH (improper) (same as check4)
    '1,6,6,8' : '28',  # CA-HN-NH-HN (improper)
    '1,1,4,10' : '29',  # CA-CA-LC-O
    '1,1,10,11' : '30',  # CA-CA-LC-LN (same as 4,6)
    # '1,6,10,11' : '31',  # CA-LC-LN-HN (same as check5)
    '1,4,10,11' : '32',  # O-LC-LN-CA (same as 5)
    '4,6,10,11' : '33',  # O-LC-LN-HN (same as 7)
    '1,1,6,11' : '34', # HN-LN-CA-CA
    '1,7,9,12' : '35', # CA-CT-OH-HO
    '4,7,9,12' : '36', # O-CT-OH-HO
    '1,1,5,6' : '37', # CA-CA-N-HN
    '1,1,7,12' : '38', # CA-CA-CT-OH
    '1,1,4,12' : '39', # CA-CA-CT-O
    '1,1,1,12' : 'check6', # CA-CA-CA-CT (same as 42)
    '1,1,2,12' : '41', # CT-CA-CA-HA
    # '1,1,1,12' : '42'  # CA-CA-CA-CT (improper) (same as check6)
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
        a_type = int(l[2])
        charge = l[3]
        x = l[4]
        y = l[5]
        z = l[6]

        atoms[atom_id] = a_type
        out.write(line)

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

        types = [a1_type,a2_type]
        types.sort()

        bond_pair = ''
        for i,t in enumerate(types):
            bond_pair += str(t)
            if i < len(types)-1:
                bond_pair += ','

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

        ang_id = l[0]
        old_type = l[1]
        a1_id = l[2]
        a2_id = l[3]
        a3_id = l[4]
        a1_type = atoms[a1_id]
        a2_type = atoms[a2_id]
        a3_type = atoms[a3_id]

        types = [a1_type,a2_type,a3_type]
        types.sort()

        angle_group = ''
        for i,t in enumerate(types):
            angle_group += str(t)
            if i < len(types)-1:
                angle_group += ','

        new_type = angles[angle_group]

        new_line = "%s %s %s %s %s\n" %(ang_id, new_type, a1_id, a2_id, a3_id)
        out.write(new_line)


    else: # write blank lines
        out.write(line)

# Dihedrals section
for line in f:

    if len(line.split()) > 0: # if line is not blank

        l = line.split()

        dih_id = l[0]
        old_type = l[1]
        a1_id = l[2]
        a2_id = l[3]
        a3_id = l[4]
        a4_id = l[5]
        a1_type = atoms[a1_id]
        a2_type = atoms[a2_id]
        a3_type = atoms[a3_id]
        a4_type = atoms[a4_id]

        types = [a1_type,a2_type,a3_type,a4_type]
        types.sort()

        dih_group = ''
        for i,t in enumerate(types):
            dih_group += str(t)
            if i < len(types)-1:
                dih_group += ','

        new_type = dihedrals[dih_group]

        new_line = "%s %s %s %s %s %s\n" %(dih_id, new_type, a1_id, a2_id, a3_id, a4_id)
        out.write(new_line)

    else: # write blank lines
        out.write(line)
