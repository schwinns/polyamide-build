#!/usr/bin/env python

# Remove uncrosslinked oligomers and fix lmps file to have proper atom numbers after removing uncrosslinked oligomers

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to remove atoms')
args = parser.parse_args()

f = open(args.lmps,'r')

# Header and Coeffs
dims = []
for line in f:

    if line.startswith('Atoms'): # go to Atoms section
        break

    elif len(line.split()) == 4: # box dimensions
        
        lo = float(line.split()[0])
        hi = float(line.split()[1])

        dims.append([lo,hi])

# Atoms section
min_z = dims[2][0]
max_z = dims[2][1]
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        a_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        charge = l[3]
        x = float(l[4])
        y = float(l[5])
        z = float(l[6])

        if z < min_z:
            min_z = z
        if z > max_z:
            max_z = z

print('Box dimensions should be:')
print('%f %f xlo xhi' %(dims[0][0], dims[0][1]))
print('%f %f ylo yhi' %(dims[1][0], dims[1][1]))
print('%f %f zlo zhi\n' %(min_z - 5, max_z + 5))