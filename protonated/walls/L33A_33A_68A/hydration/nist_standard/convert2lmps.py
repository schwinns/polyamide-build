
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f','--file',
                    help='txt file to convert to lammps format (from nist txt format)')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()


#######################################################################################
################################## GET DATA FROM TXT FILE #############################
#######################################################################################

f = open(args.file, 'r')
out = open(args.output, 'w')

# Get header info from txt file
box = f.readline()
box_x = float(box.split()[0])
box_y = float(box.split()[1])
box_z = float(box.split()[2])

n_molecules = int(f.readline().split()[0])
n_atoms = n_molecules * 3
n_bonds = n_molecules * 2
n_angles = n_molecules

# Get coordinate data, atom types from txt file
atoms = {}
for line in f:

    l = line.split()
    a_id = int(l[0])
    x = float(l[1])
    y = float(l[2])
    z = float(l[3])
    a_type = l[4]

    atoms[a_id] = {
        'coords' : [x,y,z]
    }

    if a_type == 'O':
        atoms[a_id]['type'] = 2 # OW in params.lmps
        atoms[a_id]['charge'] = -0.834
        atoms[a_id]['bonded'] = [a_id+1, a_id+2]
    elif a_type == 'H':
        atoms[a_id]['type'] = 1 # HW in params.lmps
        atoms[a_id]['charge'] = 0.417

#######################################################################################
############################## WRITE A NEW LAMMPS FILE ################################
#######################################################################################

# Header and params
header = 'Reference SPC/E Water: %d molecules\n\n' %(n_molecules)
out.write(header)

out.write('%s atoms\n' %(n_atoms))
out.write('%s bonds\n' %(n_bonds))
out.write('%s angles\n' %(n_angles))
out.write('0 dihedrals\n')
out.write('0 impropers\n\n')

# out.write('16 atom types\n')
# out.write('19 bond types\n')
# out.write('24 angle types\n')
# out.write('42 dihedral types\n\n')

out.write('2 atom types\n')
out.write('1 bond types\n')
out.write('1 angle types\n\n')

xdim = '%.6f %.6f xlo xhi\n' %(-box_x/2, box_x/2)
ydim = '%.6f %.6f ylo yhi\n' %(-box_y/2, box_y/2)
zdim = '%.6f %.6f zlo zhi\n\n' %(-box_z/2, box_z/2)

out.write(xdim + ydim + zdim)
print('Add masses and pair coefficients manually')

# Atoms section
out.write('\nAtoms\n\n')

mol_id = 0
i = 0
for atom in atoms:

    a_id = atom
    a_type = atoms[atom]['type']
    charge = atoms[atom]['charge']
    xyz = atoms[atom]['coords']
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    if i % 3 == 0:
        mol_id += 1
    
    line = ' %s %s %s %s %s %s %s\n' %(a_id,mol_id,a_type,charge,x,y,z)
    out.write(line)
    i += 1


# Bonds section
out.write('\nBonds\n\n')

b_id = 0
for atom in atoms:

    if atoms[atom]['type'] == 2:

        for a in atoms[atom]['bonded']:
            b_id += 1
            b_type = 1
            a1 = atom
            a2 = a

            line = ' %s %s %s %s\n' %(b_id, b_type, a1, a2)
            out.write(line)


# Angles section
out.write('\nAngles\n\n')

ang_id = 0
for atom in atoms:

    if atoms[atom]['type'] == 2:
        ang_id += 1
        ang_type = 1
        a1 = atom
        a2 = atoms[atom]['bonded'][0]
        a3 = atoms[atom]['bonded'][1]

        line = ' %s %s %s %s %s\n' %(ang_id, ang_type, a1, a2, a3)
        out.write(line)