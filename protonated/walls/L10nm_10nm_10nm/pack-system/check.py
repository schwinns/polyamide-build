#!/usr/bin/env python

# Fix lmps file to have proper atom numbers for packing
# Must manually change atom types in lmps Atoms section

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='pdb file to rename atoms')
args = parser.parse_args()

f = open(args.lmps,'r')

# Header and Coeffs
for line in f:

    if line.startswith('Atoms'): # go to Atoms section
        break

quads = np.zeros(8)
# Atoms section
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        atom_id = l[0]
        mol_id = l[1]
        a_type = int(l[2])
        charge = l[3]
        x = float(l[4])
        y = float(l[5])
        z = float(l[6])

        if x < 100 and y < 100 and z < 100:   # - - - 
            quads[0] += 1
        elif x > 100 and y < 100 and z < 100: # + - -
            quads[1] += 1
        elif x < 100 and y > 100 and z < 100: # - + -
            quads[2] += 1
        elif x > 100 and y > 100 and z < 100: # + + -
            quads[3] += 1
        elif x < 100 and y < 100 and z > 100: # - - +
            quads[4] += 1
        elif x < 100 and y > 100 and z > 100: # - + +
            quads[5] += 1
        elif x > 100 and y < 100 and z > 100: # + - +
            quads[6] += 1
        elif x > 100 and y > 100 and z > 100: # + + +
            quads[7] += 1


print(quads)
