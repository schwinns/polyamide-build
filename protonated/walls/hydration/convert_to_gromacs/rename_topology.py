#!/usr/bin/env python

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-t','--top',
                    help='gromacs topology file to rename')
parser.add_argument('-o','--output',default='output.top',
                    help='output top filename')
args = parser.parse_args()

f = open(args.top,'r')
out = open(args.output, 'w')

# Dictionary with common atoms --> at_nums[element_number] = element
at_nums = {
    '1' : 'H',
    '2' : 'He',
    '5' : 'B',
    '6' : 'C',
    '7' : 'N',
    '8' : 'O',
    '9' : 'F',
    '11': 'Na',
    '17': 'Cl'
}

# Skip to [ atomtypes ]
for line in f:

    if line.startswith('[ atomtypes ]'):
        out.write(line)
        break

    else:
        out.write(line)

# Build dictionaries based on [ atomtypes ]
lmps2gro = {} # lmps2gro[lmps_atom_name] = GAFF Gromacs atom type
for line in f:

    if line.startswith('[ moleculetype ]'):
        out.write(line)
        break
    
    elif line.startswith(';'): # write comments
        out.write(line)

    elif len(line.split()) > 0:

        if line.split(';') == 1:
            print('Please add GAFF atom names to [ atomtypes ] as a comment at the end of the line...')
            exit()

        gro_params = line.split('; ')[1]
        gro_name = gro_params.split()[0]
        at_num = gro_params.split()[1]

        l = line.split()
        a_type = l[0]
        b_type = l[1]
        mass = l[3]
        charge = l[4]
        ptype = l[5]
        sig = l[6]
        eps = l[7].split(';')[0]

        lmps2gro[a_type] = gro_params.split()

        new_line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t\t%s\n' %(gro_name, gro_name, at_num, mass, charge, ptype, sig, eps)
        out.write(new_line)

    else: # blank lines
        out.write(line)


# [ moleculetype ] directive
mol_names = []
for line in f:

    if line.startswith('[ atoms ]'):
        out.write(line)
        break
    
    elif line.startswith(';'): # write comments
        out.write(line)

    elif len(line.split()) > 0:

        if line.split(';') == 1:
            print('Please add molecule name to [ moleculetype ] as a comment at the end of the line...')
            exit()

        gro_params = line.split('; ')[1]
        mol_name = gro_params.split()[0]
        mol_names.append(mol_name)
        nexcl = line.split()[1].split(';')[0]

        new_line = '%s\t\t%s\n' %(mol_name, nexcl)
        out.write(new_line)

    else: # blank lines
        out.write(line)


elements = {} # elements[element] = number of element in system
for a in at_nums: # initialize elements dictionary
    element = at_nums[a]
    elements[element] = 0

# [ atoms ] directive
for line in f:

    if line.startswith('[ pairs ]'):
        out.write(line)
        break
    
    elif line.startswith(';'): # write comments
        out.write(line)

    elif len(line.split()) > 0:

        l = line.split()
        num = int(l[0])
        a_type = l[1]
        resnum = l[2]
        resname = l[3]
        cgnr = l[5]
        charge = l[6]
        mass = l[7]

        gro_type = lmps2gro[a_type][0]
        at_num = lmps2gro[a_type][1]
        element = at_nums[at_num]

        if elements[element] == 0:
            a_name = element
        else:
            a_name = element + str(elements[element])
        elements[element] += 1

        new_line = '%5i %s\t\t\t%s %s%6s\t\t%s\t%s\t%s\n' %(num, gro_type, resnum, resname, a_name, cgnr, charge, mass)
        out.write(new_line)

    else: # blank lines
        out.write(line)


# [ pairs ], [ bonds ], and [ angles ] directives --> should be unchanged
for line in f:

    if line.startswith('[ dihedrals ]'):
        out.write(line)
        break

    else:
        out.write(line)


# [ dihedrals ] directive
for line in f:

    if line.startswith('[ system ]'):
        out.write(line)
        break

    elif line.startswith(';'): # comments
        out.write(line)

    elif len(line.split()) > 8: # correct dihedrals
        out.write(line)

    elif len(line.split()) > 0: # manually changed 7 or 36 dihedrals

        l = line.split()
        a1 = int(l[0])
        a2 = int(l[1])
        a3 = int(l[2])
        a4 = int(l[3])

        new_line1 = '%7i%8i%8i%8i\t 1\t 0.0000000\t\t7.9496000\t\t1\n' %(a1,a2,a3,a4)
        new_line2 = '%7i%8i%8i%8i\t 1\t 180.0000771\t\t9.6232000\t\t2\n' %(a1,a2,a3,a4)
        out.write(new_line1)
        out.write(new_line2)

    else:
        out.write(line)


# [ system ] directive
for line in f:

    if line.startswith('[ molecules ]'):
        out.write(line)
        break

    else:
        out.write(line)


# [ molecules ] directive
n = 0
for line in f:

    if line.startswith(';'): # comments
        out.write(line)

    elif len(line.split()) > 0:

        nmols = line.split()[1]
        comp = mol_names[n]
        n += 1

        new_line = '%s\t\t%s\n' %(comp, nmols)
        out.write(new_line)

    else: # write all other lines
        out.write(line)
