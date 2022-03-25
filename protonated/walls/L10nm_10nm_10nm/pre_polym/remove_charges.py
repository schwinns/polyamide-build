#!/usr/bin/env python

# Fix lmps file to have no partial charges

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps data file to remove partial charges')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

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
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        atom_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        old_charge = l[3]
        x = l[4]
        y = l[5]
        z = l[6]

        new_charge = 0

        new_line = "%s %s %s %d %s %s %s\n" %(atom_id, mol_id, a_type, new_charge, x, y, z)
        out.write(new_line)

    else: # write blank lines
        out.write(line)


# Bonds, angle, dihedral, etc
for line in f:

    out.write(line)
