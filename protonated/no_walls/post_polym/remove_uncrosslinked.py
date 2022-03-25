#!/usr/bin/env python

# Remove uncrosslinked oligomers and fix lmps file to have proper atom numbers after removing uncrosslinked oligomers

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to remove atoms')
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
i = 1
unxlink_N = 0; unxlink_C = 0
xlink_N = 0; xlink_C = 0
term_N = 0; term_C = 0
id_key = {}
delete_id = []
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        old_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        charge = l[3]
        x = l[4]
        y = l[5]
        z = l[6]

        if mol_id == str(1): # only write the first molecule (the crosslinked molecule)
            new_line = " %d %s %s %s %s %s %s\n" %(i,mol_id,a_type,charge,x,y,z)
            out.write(new_line)

            id_key[old_id] = i # Save the conversion between old_id and new_id
            i += 1

            if a_type == str(5):
                unxlink_N += 1
            elif a_type == str(11):
                xlink_N += 1
            elif a_type == str(8):
                term_N += 1
            elif a_type == str(3):
                unxlink_C += 1
            elif a_type == str(10):
                xlink_C += 1
            elif a_type == str(12):
                term_C += 1
    
        else:
            delete_id.append(old_id) # Save the deleted atoms for later sections
            print('Molecule %s, atom type %s' %(mol_id, a_type))

    else: # write blank lines
        out.write(line)

print('%d amine groups are crosslinked, %d are uncrosslinked, and %d are terminated' %(xlink_N, unxlink_N, term_N))
print('%d carboxyl groups are crosslinked, %d are uncrosslinked, and %d are terminated' %(xlink_C, unxlink_C, term_C))
print('Degree of crosslinking is %.3f\n' %(xlink_N / (xlink_N + unxlink_N + term_N)))

print('Removing %d atoms from %s' %(len(delete_id), args.lmps))
print('%s has %d atoms' %(args.output, i-1))

# Bonds section
i = 1
for line in f:

    if line.startswith('Angles'): # go to Angles section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        b_type = l[1]
        a1_old = l[2]
        a2_old = l[3]

        if a1_old not in delete_id and a2_old not in delete_id: # if the atom is not removed (aka if it is in crosslinked molecule)

            a1_new = id_key[a1_old]
            a2_new = id_key[a2_old]

            new_line = " %d %s %s %s\n" %(i,b_type,a1_new,a2_new)
            out.write(new_line)
            i += 1

    else: # write blank lines
        out.write(line)

print('%s has %d bonds' %(args.output, i-1))

# Angles section
i = 1
for line in f:

    if line.startswith('Dihedrals'): # go to Dihedrals section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        ang_type = l[1]
        a1_old = l[2]
        a2_old = l[3]
        a3_old = l[4]

        if a1_old not in delete_id and a2_old not in delete_id and a3_old not in delete_id:

            a1_new = id_key[a1_old]
            a2_new = id_key[a2_old]
            a3_new = id_key[a3_old]

            new_line = " %d %s %s %s %s\n" %(i,ang_type,a1_new,a2_new,a3_new)
            out.write(new_line)
            i += 1

    else: # write blank lines
        out.write(line)

print('%s has %d angles' %(args.output, i-1))

# Dihedrals section
i = 1
for line in f:

    if line.startswith('Impropers'): # go to Impropers section
        out.write(line)
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        d_type = l[1]
        a1_old = l[2]
        a2_old = l[3]
        a3_old = l[4]
        a4_old = l[5]

        if a1_old not in delete_id and a2_old not in delete_id and a3_old not in delete_id and a4_old not in delete_id:

            a1_new = id_key[a1_old]
            a2_new = id_key[a2_old]
            a3_new = id_key[a3_old]
            a4_new = id_key[a4_old]

            new_line = " %d %s %s %s %s %s\n" %(i,d_type,a1_new,a2_new,a3_new,a4_new)
            out.write(new_line)
            i += 1

    else: # write blank lines
        out.write(line)

print('%s has %d dihedrals' %(args.output, i-1))

# Impropers section
i = 1
for line in f:

    if len(line.split()) > 0: # if line is not blank

        l = line.split()

        i_type = l[1]
        a1_old = l[2]
        a2_old = l[3]
        a3_old = l[4]
        a4_old = l[5]

        if a1_old not in delete_id and a2_old not in delete_id and a3_old not in delete_id and a4_old not in delete_id:

            a1_new = id_key[a1_old]
            a2_new = id_key[a2_old]
            a3_new = id_key[a3_old]
            a4_new = id_key[a4_old]

            new_line = " %d %s %s %s %s %s\n" %(i,i_type,a1_new,a2_new,a3_new,a4_new)
            out.write(new_line)
            i += 1

    else: # write blank lines
        out.write(line)

print('%s has %d impropers' %(args.output, i-1))