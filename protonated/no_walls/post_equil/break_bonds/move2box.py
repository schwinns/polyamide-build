# Fix lmps file to place all atoms within the periodic box

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to move atoms')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

f = open(args.lmps, 'r')
out = open(args.output, 'w')

# Get box lengths and write everything before atoms
for line in f:

    if line.startswith('Atoms'):
        out.write(line)
        break

    elif len(line.split()) == 4:

        l = line.split()
        lo = float(l[0])
        hi = float(l[1])

        box = hi - lo
        out.write(line)

    else:
        out.write(line)

print('Box dimensions (Angstroms):')
print('\t%f %f' %(lo,hi))
print('Box length:')
print('\t%f\n' %(hi-lo))

# Move all atoms within periodic box
for line in f:

    if line.startswith('Bonds'): # end of atoms section
        out.write(line)
        break

    elif len(line.split()) > 0: # if not blank

        l = line.split()
        a_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        charge = l[3]
        x = float(l[4])
        y = float(l[5])
        z = float(l[6])

        dims = [x,y,z]
        for i,d in enumerate(dims):
            if d > hi:
                dims[i] = d - box
            elif d < lo:
                dims[i] = d + box

        new_line = '%s %s %s %s %f %f %f\n' %(a_id, mol_id, a_type, charge, dims[0], dims[1], dims[2])
        out.write(new_line)

    else: # write blank lines
        out.write(line)

# Write all other information
for line in f:
    out.write(line)