#!/usr/bin/env python

# Combine lmps files with water and PA membrane

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-pa','--lmps_PA',
                    help='lammps file with PA membrane')
parser.add_argument('-w','--lmps_water',
                    help='lmps file with water box')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

###############################################################################
############################ SAVE WATER TOPOLOGY ##############################
###############################################################################

# lmp = open(args.lmps_water, 'r')
lmp = open('final.lmps', 'r')

# Read in all information from water lammps file to add to end of PA membrane file
masses = {}
pair_coeffs = {}
for line in lmp:

    if line.startswith('Atoms'):
        break

    elif len(line.split()) == 2:
        # Number of atoms, bonds, etc
        if 'atoms' in line.split():
            n_atoms = int(line.split()[0])
        elif 'bonds' in line.split():
            n_bonds = int(line.split()[0])
        elif 'angles' in line.split():
            n_angles = int(line.split()[0])
        elif 'dihedrals' in line.split():
            n_dihedrals = int(line.split()[0])
        
        # Masses
        else:
            masses[line.split()[0]] = line.split()[1]

    elif len(line.split()) == 3:
        # Types of atoms, bonds, etc
        if 'atom' in line.split():
            a_types = int(line.split()[0])
        elif 'bond' in line.split():
            b_types = int(line.split()[0])
        elif 'angle' in line.split():
            ang_types = int(line.split()[0])
        elif 'dihedral' in line.split():
            dih_types = int(line.split()[0])

        # Pair Coeffs
        else:
            pair_coeffs[line.split()[0]] = line.split()[1:]

    elif len(line.split()) == 4 and '#' not in line.split(): # box dimensions
        if 'xlo' in line.split():
            xlo = float(line.split()[0])
            xhi = float(line.split()[1])
        elif 'ylo' in line.split():
            ylo = float(line.split()[0])
            yhi = float(line.split()[1])
        elif 'zlo' in line.split():
            zlo = float(line.split()[0])
            zhi = float(line.split()[1])

# Build dictionary with atomistic topology -- atoms section
atoms = {}
for line in lmp:
    if line.startswith('Bonds'):
        vel = False
        break

    elif line.startswith('Velocities'):
        vel = True
        break

    elif len(line.split()) > 0:
        
        l = line.split()
        a_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        charge = l[3]
        x = float(l[4])
        y = float(l[5])
        z = float(l[6])

        atoms[a_id] = {
            'mol'    : mol_id,
            'type'   : a_type,
            'charge' : charge,
            'coords' : np.array([x,y,z]),
            'bonded' : [],
            'delete' : False
        }

# Skip to bonding section
if vel:
    for line in lmp:
        if line.startswith('Bonds'):
            break

# Build dictionary with atomistic topology -- bonds section
bonds = {}
for line in lmp:

    if line.startswith('Angles'):
        break

    elif len(line.split()) > 0:

        l = line.split()
        b_id = l[0]
        b_type = l[1]
        a1 = l[2]
        a2 = l[3]

        atoms[a1]['bonded'].append(a2)
        atoms[a2]['bonded'].append(a1)

        bonds[b_id] = {
            'type'  : b_type,
            'atoms' : [a1,a2]
        }

# Build dictionary with atomistic topology -- angles section
angles = {}
for line in lmp:

    if line.startswith('Dihedrals'):
        break

    elif len(line.split()) > 0:

        l = line.split()
        ang_id = l[0]
        ang_type = l[1]
        a1 = l[2]
        a2 = l[3]
        a3 = l[4]

        angles[ang_id] = {
            'type'  : ang_type,
            'atoms' : [a1,a2,a3]
        }

# Build dictionary with atomistic topology -- dihedrals section
dihedrals = {}
for line in lmp:

    if line.startswith('Impropers'):
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
            'type'  : dih_type,
            'atoms' : [a1,a2,a3,a4]
        }

if len(atoms) != n_atoms:
    print('Numbers of atoms do not match...')
    exit()
if len(bonds) != n_bonds:
    print('Numbers of bonds do not match...')
    exit()
if len(angles) != n_angles:
    print('Numbers of angles do not match...')
    exit()
if len(dihedrals) != n_dihedrals:
    print('Numbers of dihedrals do not match...')
    exit()

###############################################################################
######################### ADD WATERS TO PA MEMBRANE ###########################
###############################################################################

# out = open(args.output, 'w')
# lmp = open(args.lmps_PA, 'r')
out = open('output.lmps' ,'w')
lmp = open('post_equil.lmps', 'r')

# Write header
m = 0
for line in lmp:

    if line.startswith('Masses'):
        out.write(line + '\n')
        break
    
    elif len(line.split()) == 2:
        # Number of atoms, bonds, etc
        if 'atoms' in line.split():
            n = 2*n_atoms + int(line.split()[0])
            n_PA = int(line.split()[0])
            out.write('%d atoms\n' %(n))
        elif 'bonds' in line.split():
            n = 2*n_bonds + int(line.split()[0])
            out.write('%d bonds\n' %(n))
        elif 'angles' in line.split():
            n = 2*n_angles + int(line.split()[0])
            out.write('%d angles\n' %(n))
        elif 'dihedrals' in line.split():
            n = 2*n_dihedrals + int(line.split()[0])
            out.write('%d dihedrals\n' %(n))

    elif len(line.split()) == 3:
        # Types of atoms, bonds, etc
        if 'atom' in line.split():
            out.write('%s atom types\n' %(a_types) )
        elif 'bond' in line.split():
            out.write('%s bond types\n' %(b_types) )
        elif 'angle' in line.split():
            out.write('%s angle types\n' %(ang_types) )
        elif 'dihedral' in line.split():
            out.write('%s dihedral types\n' %(dih_types) )

    elif len(line.split()) == 4 and '#' not in line.split(): # box dimensions
        if 'xlo' in line.split():
            PA_xlo = float(line.split()[0])
            PA_xhi = float(line.split()[1])
            out.write(line)
        elif 'ylo' in line.split():
            PA_ylo = float(line.split()[0])
            PA_yhi = float(line.split()[1])
            out.write(line)
        elif 'zlo' in line.split(): # 4 nm above and below PA membrane
            PA_zlo = float(line.split()[0])
            PA_zhi = float(line.split()[1])

            out.write('%f %f zlo zhi\n' %(PA_zlo - (zhi-zlo), PA_zhi + (zhi-zlo)) ) # for now, just a single water box on each side
            #out.write('%f %f zlo zhi\n' %(PA_zlo - 40, PA_zhi + 40) )

    else:
        out.write(line)

print('PA box dimensions:')
print('\t%.4f %.4f xlo xhi\n\t%.4f %.4f ylo yhi\n\t%.4f %.4f zlo zhi' %(PA_xlo, PA_xhi, PA_ylo, PA_yhi, PA_zlo, PA_zhi) )

for mass in masses: # write masses from water file
    out.write('%s %s\n' %(mass, masses[mass]) )

out.write('\nPair Coeffs\n\n')
for p in pair_coeffs: # write pair coeffs from water file
    params = pair_coeffs[p]
    out.write('%s %s %s\n' %(p, params[0], params[1]) )
    

# Write Atoms section and save atom ids for renumbering atoms
out.write('\n')
atom_ids = {} # atom_ids[old_id] = new_id
i = 1
for line in lmp: # skip to atoms section

    if line.startswith('Atoms'):
        out.write(line + '\n')
        break

for line in lmp: 

    if line.startswith('Bonds'):
        vel = False
        break

    elif line.startswith('Velocities'):
        vel = True
        break

    elif len(line.split()) > 0: # write all the PA atoms unchanged
        
        l = line.split()
        old_id = l[0]
        mol_id = l[1]
        a_type = l[2]
        charge = l[3]
        x = l[4]
        y = l[5]
        z = l[6]


        new_line = '%d %s %s %s %s %s %s\n' %(i, mol_id, a_type, charge, x, y, z)
        out.write(new_line)

        atom_ids[old_id] = i
        i += 1


# TOP RESERVOIR pt 1 (from PA to full water box in z)
n_top = 0
n_waters = 0
for atom in atoms:

    a_id = n_PA + int(atom)
    mol_id = int(atoms[atom]['mol'])
    mol_id += 1 # add the PA membrane
    a_type = atoms[atom]['type']
    charge = atoms[atom]['charge']
    xyz = atoms[atom]['coords']

    # shift all coordinates to edge of PA cross-section
    x_shift = xyz[0] - xlo + PA_xlo
    y_shift = xyz[1] - ylo + PA_ylo

    # shift all coordinates above PA membrane
    z_shift = xyz[2] - zlo + PA_zhi

    # if outside PA cross section
    if x_shift > PA_xhi:
        atoms[atom]['delete'] = True
    elif y_shift > PA_yhi:
        atoms[atom]['delete'] = True

    # if atoms[atom]['delete']:
    #     for a in atoms[atom]['bonded']: # delete bonded atoms
    #         atoms[a]['delete'] = True

    if a_type == '15' and not atoms[atom]['delete']: # save waters where oxygen is within PA membrane cross section
        n_top += 1
        new_line = '%d %s %s %s %s %s %s\n' %(i, mol_id, a_type, charge, x_shift, y_shift, z_shift)
        out.write(new_line)
        atom_ids[a_id] = i
        i += 1

        for a in atoms[atom]['bonded']: # save hydrogens bonded to good oxygens
            


print('i',i)
print('PA atoms',n_PA)
print('Top reservoir atoms',n_top)
print('Top waters',n_waters)
# # TOP RESERVOIR pt 2 (from water box to +40 A in z)
# top_pt1_atom = n_top
# top_pt1_mol = n_top / 3
# for atom in atoms:

#     a_id = n_PA + top_pt1_atom + int(atom)
#     mol_id = int(atoms[atom]['mol']) + top_pt1_mol
#     mol_id += 1 # add the PA membrane
#     a_type = atoms[atom]['type']
#     charge = atoms[atom]['charge']
#     xyz = atoms[atom]['coords']

#     # shift all coordinates to edge of PA cross-section
#     x_shift = xyz[0] - xlo + PA_xlo
#     y_shift = xyz[1] - ylo + PA_ylo

#     # shift all coordinates above PA membrane + water box
#     z_shift = xyz[2] - zlo + PA_zhi + (zhi - zlo)

#     # if outside periodic box
#     if z_shift > PA_zhi + 40:
#         atoms[atom]['delete'] = True

#     if atoms[atom]['delete']: 
#         for a in atoms[atom]['bonded']: # delete bonded atoms
#             atoms[a]['delete'] = True

#     if not atoms[atom]['delete']: # carry over the deleted atoms from before
#         new_line = '%d %s %s %s %s %s %s\n' %(i, mol_id, a_type, charge, x_shift, y_shift, z_shift)
#         out.write(new_line)
#         atom_ids[a_id] = i
#         i += 1
#         n_top += 1

# BOTTOM RESERVOIR pt 1 (from PA to -full water box in z)
for atom in atoms:

    a_id = n_PA + int(atom) + n_top
    mol_id = int(atoms[atom]['mol']) + n_waters
    mol_id += 1 # add the PA membrane
    a_type = atoms[atom]['type']
    charge = atoms[atom]['charge']
    xyz = atoms[atom]['coords']

    # shift all coordinates to edge of PA cross-section
    x_shift = xyz[0] - xlo + PA_xlo
    y_shift = xyz[1] - ylo + PA_ylo

    # shift all coordinates below PA membrane
    z_shift = xyz[2] - zlo - PA_zlo - (zhi-zlo)

    if not atoms[atom]['delete']:
        new_line = '%d %s %s %s %s %s %s\n' %(i, mol_id, a_type, charge, x_shift, y_shift, z_shift)
        out.write(new_line)
        atom_ids[a_id] = i
        i += 1
        n_top += 1

# # BOTTOM RESERVOIR pt 2 (from water box to -40 A in z)
# top_pt1_atom = n_top
# top_pt1_mol = n_top / 3
# for atom in atoms:

#     a_id = n_PA + top_pt1_atom + int(atom)
#     mol_id = int(atoms[atom]['mol']) + top_pt1_mol
#     mol_id += 1 # add the PA membrane
#     a_type = atoms[atom]['type']
#     charge = atoms[atom]['charge']
#     xyz = atoms[atom]['coords']

#     # shift all coordinates to edge of PA cross-section
#     x_shift = xyz[0] - xlo + PA_xlo
#     y_shift = xyz[1] - ylo + PA_ylo

#     # shift all coordinates above PA membrane + water box
#     z_shift = xyz[2] - zlo + PA_zhi + (zhi - zlo)

#     # if outside periodic box
#     if z_shift > PA_zhi + 40:
#         atoms[atom]['delete'] = True

#     if atoms[atom]['delete']: 
#         for a in atoms[atom]['bonded']: # delete bonded atoms
#             atoms[a]['delete'] = True

#     if not atoms[atom]['delete']: # carry over the deleted atoms from before
#         new_line = '%d %s %s %s %s %s %s\n' %(i, mol_id, a_type, charge, x_shift, y_shift, z_shift)
#         out.write(new_line)
#         atom_ids[a_id] = i
#         i += 1
#         n_top += 1

n_atoms = i - 1

# skip velocities if present
if vel:
    for line in lmp: 
        if line.startswith('Bonds'):
            out.write(line)
            break

# Write Bonds section
b = 1
for line in lmp: # write the PA bonds with new atom numbers

    if line.startswith('Angles'):
        break

    elif len(line.split()) > 0:
        
        l = line.split()
        b_id = l[0]
        b_type = l[1]
        a1_old = l[2]
        a2_old = l[3]

        a1_new = atom_ids[a1_old]
        a2_new = atom_ids[a2_old]

        new_line = '%s %s %s %s\n' %(b, b_type, a1_new, a2_new)
        out.write(new_line)

for bond in bonds: # write water bonds with new atom numbers for TOP RESERVOIR

    a1_old = bonds[bond]['atoms'][0] 
    a2_old = bonds[bond]['atoms'][1]

    if not atoms[a1_old]['delete'] and not atoms[a2_old]['delete']:

        a1_new = atom_ids[a1_old]
        a2_new = atom_ids[a2_old]




# Add new bonds for waters
for mol in waters:

    OW = waters[mol][0]

    for a in range(2):

        n_bonds += 1
        new_line = "%d %d %s %s\n" %(n_bonds, 19, OW, waters[mol][a+1])
        out.write(new_line)

out.write('\nAngles\n\n')

# Write Angles section
for line in lmp:

    if line.startswith('Dihedrals'):
        break

    elif len(line.split()) > 0:
        n_angles = int(line.split()[0])
        out.write(line)


# Add new angles for waters
for mol in waters:

    n_angles += 1
    new_line = "%d %d %s %s %s\n" %(n_angles, 24, waters[mol][0], waters[mol][1], waters[mol][2])
    out.write(new_line)

out.write('\nDihedrals\n')

# Write Dihedrals section
for line in lmp:
    out.write(line)

print('\n%s atoms' %(n_atoms))
print('%d bonds' %(n_bonds))
print('%d angles' %(n_angles))