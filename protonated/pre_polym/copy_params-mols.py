#!/usr/bin/env python

# Fix lmps file to have parameters and bonding information

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-p','--params',
                    help='lmps file with parameters')
parser.add_argument('-d','--data',
                    help='lmps file with coordinate data')
parser.add_argument('-b','--bonds',
                    help='lmps file with the bonding information')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

params = open(args.params,'r')
out = open(args.output, 'w')

# Header and Coeffs
for line in params:

    if line.startswith('Atoms'): # go to Atoms section
        out.write(line)
        break

    else: # everything before atoms
        out.write(line)

out.write('\n')
params.close()
data = open(args.data, 'r')

# Atoms section
check = False
ts = 0
for line in data: # find the last timestep (assuming 0 and last are saved)

    if line.startswith('ITEM: TIMESTEP'):
        check = True
    
    elif check:
        ts = int(line.split()[0])
        check = False

    if ts > 0:
        break

for line in data: # skip info before atom coordinates

    if line.startswith('ITEM: ATOMS'): # go to Bonds section
        break

for line in data: # write all the atom coordinates

    out.write(line)

out.write('\n')
data.close()
bonds = open(args.bonds, 'r')

# Bonds, angle, dihedral information
for line in bonds: # skip everything before the bonding info

    if line.startswith('Bonds'):
        out.write(line)
        break

for line in bonds: # write all the bond information

    out.write(line)

