# Fix lmps file to break bonds across periodic boundaries

import argparse
import numpy as np
import random

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to move atoms')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

f = open(args.lmps, 'r')
out = open(args.output, 'w')

#######################################################################################
############################ INITIAL INFORMATION FROM HEADER ##########################
#######################################################################################

# Get box lengths and write masses and pair coeffs
for line in f:

    if line.startswith('Atoms'):
        out.write(line)
        out.write('\n')
        break

    elif line.startswith('Bond Coeffs'):
        break

    elif len(line.split()) == 4:

        l = line.split()
        lo = float(l[0])
        hi = float(l[1])

        box = hi - lo
        out.write(line)

    else:
        if '#' in line.split():
            new_line = line.split('#')[0]
            out.write(new_line + '\n')
        else:
            out.write(line)

for line in f: # skip to atoms

    if line.startswith('Atoms'):
        out.write(line)
        out.write('\n')
        break

#######################################################################################
######################### CREATE DICTIONARIES WITH FULL TOPOLOGY ######################
#######################################################################################

# Create dictionaries with full topology
atoms = {}
unxlink_N = 0; unxlink_C = 0
xlink_N = 0; xlink_C = 0
term_N = 0; term_C = 0
for line in f: # build dictionary of all atoms

    if line.startswith('Bonds'): # end of atoms section
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

        atoms[a_id] = {
            'molecule' : mol_id,
            'type'   : a_type,
            'charge' : charge,
            'coords' : np.array([x,y,z])
        }

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

bonds = {}
for line in f: # build dictionary of bonds

    if line.startswith('Angles'): # end of bonds section
        break
    
    elif len(line.split()) > 0:

        l = line.split()
        b_id = l[0]
        b_type = l[1]
        a1 = l[2]
        a2 = l[3]
        
        bonds[b_id] = {
            'type' : b_type,
            'atoms' : [a1,a2],
        }

angles = {}
for line in f: # build angles dictionary

    if line.startswith('Dihedrals'): # end of angles section
        break
    
    elif len(line.split()) > 0:

        l = line.split()
        ang_id = l[0]
        ang_type = l[1]
        a1 = l[2]
        a2 = l[3]
        a3 = l[4]

        angles[ang_id] = {
            'type' : ang_type,
            'atoms' : [a1,a2,a3],
        }

dihedrals = {}
for line in f: # build dihedrals dictionary

    if line.startswith('Impropers'): # end of dihedrals section
        break
    
    elif len(line.split()) > 0:

        l = line.split()
        dih_id = l[0]
        dih_type = l[1]
        a1 = l[2]
        a2 = l[3]
        a3 = l[4]
        a4 = l[5]

        dihedrals[dih_id] = {
            'type' : dih_type,
            'atoms' : [a1,a2,a3,a4],
        }

#######################################################################################
######################## SHIFT ATOMS SO BONDS ARE NOT ACROSS BOX ######################
#######################################################################################

# Find all the bonds that cross the periodic boundary in the z direction
print()
broken = 1
i = 1
while broken != 0:

    # Shift atoms so bonds are not too large
    for bond in bonds:

        a1 = random.choice([0,1])
        if a1 == 0:
            a2 = 1
        else:
            a2 = 0


        a1 = bonds[bond]['atoms'][a1]
        a2 = bonds[bond]['atoms'][a2]

        a1_xyz = atoms[a1]['coords']
        a2_xyz = atoms[a2]['coords']

        distance = np.abs(a1_xyz[2] - a2_xyz[2]) # distance in z direction
        bonds[bond]['distance'] = distance

        if distance > box/2:
        
            if a1_xyz[2] > lo + box/2: # if a1 in top half of box, move a1 down
                atoms[a1]['coords'][2] -= box
            elif a1_xyz[2] < hi - box/2: # if a1 in bottom half of box, move a1 up
                atoms[a1]['coords'][2] += box

    # update the number of broken bonds
    broken = 0
    for bond in bonds:

        a1 = bonds[bond]['atoms'][0]
        a2 = bonds[bond]['atoms'][1]

        a1_xyz = atoms[a1]['coords']
        a2_xyz = atoms[a2]['coords']

        distance = np.abs(a1_xyz[2] - a2_xyz[2]) # distance in z direction
        bonds[bond]['distance'] = distance

        if distance > box/2:
            broken += 1
            
    if i % 100 == 0 or broken == 0:
        print('%d broken bonds in %d iterations' %(broken, i))
    
    i += 1
    

#######################################################################################
####################### WRITE A NEW LAMMPS FILE WITH SHIFTED BONDS ####################
#######################################################################################

# Header and coefficient were written above

maxes = [-100,-100,-100]
mins = [100,100,100]
# Atoms section
for atom in atoms:

    a_id = atom
    mol_id = atoms[atom]['molecule']
    a_type = atoms[atom]['type']
    charge = atoms[atom]['charge']
    x = atoms[atom]['coords'][0]
    y = atoms[atom]['coords'][1]
    z = atoms[atom]['coords'][2]

    xyz = [float(x),float(y),float(z)]
    for d in range(len(xyz)):
        if xyz[d] > maxes[d]:
            maxes[d] = xyz[d]
        elif xyz[d] < mins[d]:
            mins[d] = xyz[d]

    new_line = '%s %s %s %s %f %f %f\n' %(a_id, mol_id, a_type, charge, x, y, z)
    out.write(new_line)

# Bonds section
out.write('\nBonds\n\n')
b = 1
for bond in bonds: # need to renumber bonds because some are deleted

    b_type = bonds[bond]['type']
    a1 = bonds[bond]['atoms'][0]
    a2 = bonds[bond]['atoms'][1]

    new_line = '%d %s %s %s\n' %(b, b_type, a1, a2)
    out.write(new_line)
    b += 1

# Angles section
out.write('\nAngles\n\n')
ang = 1
for angle in angles: # need to renumber angles because some are deleted

    ang_type = angles[angle]['type']
    a1 = angles[angle]['atoms'][0]
    a2 = angles[angle]['atoms'][1]
    a3 = angles[angle]['atoms'][2]

    new_line = '%d %s %s %s %s\n' %(ang, ang_type, a1, a2, a3)
    out.write(new_line)
    ang += 1

# Dihedrals section
out.write('\nDihedrals\n\n')
dih = 1
for dihedral in dihedrals: # need to renumber dihedrals because some are deleted

    dih_type = dihedrals[dihedral]['type']
    a1 = dihedrals[dihedral]['atoms'][0]
    a2 = dihedrals[dihedral]['atoms'][1]
    a3 = dihedrals[dihedral]['atoms'][2]
    a4 = dihedrals[dihedral]['atoms'][3]

    new_line = '%d %s %s %s %s %s\n' %(dih, dih_type, a1, a2, a3, a4)
    out.write(new_line)
    dih += 1

print('Box dimensions should be:')
print('\t%f %f xlo xhi' %(mins[0], maxes[0]))
print('\t%f %f ylo yhi' %(mins[1], maxes[1]))
print('\t%f %f zlo zhi\n' %(mins[2], maxes[2]))

print('%d amine groups are crosslinked, %d are uncrosslinked, and %d are terminated' %(xlink_N, unxlink_N, term_N))
print('%d carboxyl groups are crosslinked, %d are uncrosslinked, and %d are terminated' %(xlink_C, unxlink_C, term_C))
print('Degree of crosslinking is %.3f\n' %(xlink_N / (xlink_N + unxlink_N + term_N)))