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

# Get box lengths and write everything before atoms
for line in f:

    if line.startswith('Atoms'):
        out.write(line)
        out.write('\n')
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

#######################################################################################
######################### CREATE DICTIONARIES WITH FULL TOPOLOGY ######################
#######################################################################################

# Create dictionaries with full topology
atoms = {}
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
            'coords' : np.array([x,y,z]),
            'bonded' : [],
            'bonds' : {}
        }

bonds = {}
for line in f: # build dictionary of bonds and add bond info to atom dictionary

    if line.startswith('Angles'): # end of bonds section
        break
    
    elif len(line.split()) > 0:

        l = line.split()
        b_id = l[0]
        b_type = l[1]
        a1 = l[2]
        a2 = l[3]

        atoms[a1]['bonded'].append(a2)
        atoms[a2]['bonded'].append(a1)

        atoms[a1]['bonds'][b_id] = a2
        atoms[a2]['bonds'][b_id] = a1
        
        bonds[b_id] = {
            'type' : b_type,
            'atoms' : [a1,a2],
            'broken' : False,
            'delete' : False
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
            'delete' : False
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
            'delete' : False
        }

#######################################################################################
############################ LOCATE BONDS CROSSING Z-BOUNDARY #########################
#######################################################################################

# Find all the bonds that cross the periodic boundary in the z direction
broken = 0
for bond in bonds:

    a1 = bonds[bond]['atoms'][0]
    a2 = bonds[bond]['atoms'][1]

    a1_xyz = atoms[a1]['coords']
    a2_xyz = atoms[a2]['coords']

    distance = np.abs(a1_xyz[2] - a2_xyz[2]) # distance in z direction
    bonds[bond]['distance'] = distance

    if distance > box/2:

        bonds[bond]['broken'] = True
        if bonds[bond]['type'] == '5': # mark the bonds to be deleted
            bonds[bond]['delete'] = True
        else:
            broken += 1

print('%d broken bonds\n' %(broken))

#######################################################################################
############################ LOCATE NEAREST LC-LN CROSSLINK BONDS #####################
#######################################################################################

for b in bonds:

    if int(b) % 100 == 0 or int(b) == len(bonds):
        print('{:.1%} complete with finding nearest crosslinked bond'.format(int(b) / len(bonds)))

    # Find nearest crosslinking group (10, 11) to the broken bond
    crosslinking = ['10','11']
    checked_atoms = []
    checked_bonds = []
    broken = bonds[b]['broken']
    uncrosslink = b
    while broken and not bonds[b]['delete']:

        # check all atoms bonded to the first atom in uncrosslinked bond
        a1 = random.choice([0,1])
        if a1 == 0:
            a2 = 1
        else:
            a2 = 0
        first_atom = bonds[uncrosslink]['atoms'][a1]
        bonded_atoms = atoms[first_atom]['bonded']
        checked_atoms.append(first_atom)
        
        for a in bonded_atoms:
            if a not in checked_atoms:
                next_atom = atoms[a]

                if next_atom['type'] in crosslinking:
                    for bid in atoms[a]['bonds']: # find the crosslinking bond
                        if bonds[bid]['type'] == '5': 
                            bonds[bid]['broken'] = True
                            bonds[bid]['delete'] = True
                            # print('Deleting bond %s between atoms %s and %s' %(bid, bonds[bid]['atoms'][0], bonds[bid]['atoms'][1]))
                            break
                    broken = False
                    break
                else:
                    checked_atoms.append(a)

        if not broken:
            break

        # check all atoms bonded to the second atom in uncrosslinked bond
        second_atom = bonds[uncrosslink]['atoms'][a2]
        bonded_atoms = atoms[second_atom]['bonded']
        
        for a in bonded_atoms:
            if a not in checked_atoms:
                next_atom = atoms[a]

                if next_atom['type'] in crosslinking: # if the next atom is a crosslinking atom
                    for bid in atoms[a]['bonds']: # find the crosslinking bond
                        if bonds[bid]['type'] == '5': 
                            bonds[bid]['broken'] = True
                            bonds[bid]['delete'] = True
                            # print('Deleting bond %s between atoms %s and %s' %(bid, bonds[bid]['atoms'][0], bonds[bid]['atoms'][1]))
                            break
                    broken = False
                    break
                else:
                    checked_atoms.append(a)

        if not broken:
            break

        checked_bonds.append(uncrosslink)
        b0 = uncrosslink
        i = 0
        while b0 in checked_bonds:
            # Randomly choose a new bond on randomly chosen atom in current bond
            a = random.choice([0,1])
            rand_atom = bonds[uncrosslink]['atoms'][a]
            
            b_options = [b for b in atoms[rand_atom]['bonds']]
            b0 = random.choice(b_options)
            if i == 50: # force break if we get stuck in loop and startover
                b0 = b
                checked_bonds = []
                checked_atoms = []
                print('Force break at iteration 50. Retrying at bond %s' %(b))
                break
            else:
                i += 1

        uncrosslink = b0
    

# Angle and dihedral types to change
to_change = {
    'angles' : {
        '15' : '3', # CA-LC-O  --> CA-C-O
        '16' : '4', # CA-CA-LC --> CA-CA-C
        '17' : '9', # LN-CA-CA --> N-CA-CA
        '18' : '10' # CA-LN-HN --> CA-N-HN
    },
    'dihedrals' : {
        '10' : '9', # CA-CA-CA-LC --> CA-CA-CA-C
        '13' : '12', # LC-CA-CA-HA --> C-CA-CA-HA
        '23' : '22', # CA-CA-CA-LC (imp) --> CA-CA-CA-C (imp)
        '29' : '3', # CA-CA-LC-O --> CA-CA-C-O
        '15' : '14', # LN-CA-CA-CA --> N-CA-CA-CA
        '17' : '16', # LN-CA-CA-HA --> N-CA-CA-HA
        '26' : '25', # CA-CA-CA-LN (imp) --> CA-CA-CA-N (imp)
        '34' : '37' # HN-LN-CA-CA --> HN-N-CA-CA
    }
}

# Mark the angles and dihedrals that should be deleted/changed in final lmps file
for bond in bonds:

    if bonds[bond]['delete']:

        a1 = bonds[bond]['atoms'][0]
        a2 = bonds[bond]['atoms'][1]

        # change LC and LN types to their uncrosslinked versions (C, N)
        if atoms[a1]['type'] == '10':
            C = a1
            N = a2
            atoms[C]['type'] = '3'
            atoms[N]['type'] = '4'
        elif atoms[a1]['type'] == '11':
            N = a1
            C = a2
            atoms[N]['type'] = '4'
            atoms[C]['type'] = '3'
       

        for angle in angles:
            if C in angles[angle]['atoms'] and N in angles[angle]['atoms']: # both in the angle, delete
                angles[angle]['delete'] = True
            elif C in angles[angle]['atoms'] or N in angles[angle]['atoms']: # if one in the angle, change
                old_type = angles[angle]['type']
                new_type = to_change['angles'][old_type]
                angles[angle]['type'] = new_type


        for dih in dihedrals:
            if C in dihedrals[dih]['atoms'] and N in dihedrals[dih]['atoms']:
                dihedrals[dih]['delete'] = True
            elif C in dihedrals[dih]['atoms'] and N in dihedrals[dih]['atoms']:
                old_type = dihedrals[dih]['type']
                new_type = to_change['dihedrals'][old_type]
                dihedrals[dih]['type'] = new_type

#######################################################################################
########### CREATE NEW ATOMS/BONDS/ANGLES/DIHEDRALS FOR TERMINATION GROUPS ############
#######################################################################################

# for bond in bonds:

#     if bonds[bond]['delete']:

#         n_atoms = len(atoms)

#         a1 = bonds[bond]['atoms'][0]
#         a2 = bonds[bond]['atoms'][1]

#         if atoms[a1]['type'] == '10':
#             LC = a1
#             LN = a2
#         else:
#             LC = a2
#             LN = a1

#         HN = atoms[LC]
#         OH = atoms[LN]

#         HN['type'] = '6'
#         old_bonded = HN['bonded']
#         HN['bonded'] = []
#         for a in old_bonded:
#             if atoms[a]['type'] == '11':
#                 HN['bonded'].append(a)
        
#         atoms[str(n_atoms + 1)] = HN

#         print('HN:',HN)
#         exit()
        

#######################################################################################
######################## WRITE A NEW LAMMPS FILE WITH BROKEN BONDS ####################
#######################################################################################

# Header and coefficient were written above

# Atoms section
a = 0
for atom in atoms:

    a_id = atom
    mol_id = atoms[atom]['molecule']
    a_type = atoms[atom]['type']
    charge = atoms[atom]['charge']
    x = atoms[atom]['coords'][0]
    y = atoms[atom]['coords'][1]
    z = atoms[atom]['coords'][2]

    new_line = '%s %s %s %s %f %f %f\n' %(a_id, mol_id, a_type, charge, x, y, z)
    out.write(new_line)
    a += 1

print('%s atoms' %(a))

# Bonds section
out.write('\nBonds\n\n')
b = 1
for bond in bonds: # need to renumber bonds because some are deleted

    if not bonds[bond]['delete']:
    
        b_type = bonds[bond]['type']
        a1 = bonds[bond]['atoms'][0]
        a2 = bonds[bond]['atoms'][1]

        new_line = '%d %s %s %s\n' %(b, b_type, a1, a2)
        out.write(new_line)
        b += 1

print('%s bonds' % (b-1))

# Angles section
out.write('\nAngles\n\n')
ang = 1
for angle in angles: # need to renumber angles because some are deleted

    if not angles[angle]['delete']:
    
        ang_type = angles[angle]['type']
        a1 = angles[angle]['atoms'][0]
        a2 = angles[angle]['atoms'][1]
        a3 = angles[angle]['atoms'][2]

        new_line = '%d %s %s %s %s\n' %(ang, ang_type, a1, a2, a3)
        out.write(new_line)
        ang += 1

print('%d angles' %(ang-1))

# Dihedrals section
out.write('\nDihedrals\n\n')
dih = 1
for dihedral in dihedrals: # need to renumber dihedrals because some are deleted

    if not dihedrals[dihedral]['delete']:
    
        dih_type = dihedrals[dihedral]['type']
        a1 = dihedrals[dihedral]['atoms'][0]
        a2 = dihedrals[dihedral]['atoms'][1]
        a3 = dihedrals[dihedral]['atoms'][2]
        a4 = dihedrals[dihedral]['atoms'][3]

        new_line = '%d %s %s %s %s %s\n' %(dih, dih_type, a1, a2, a3, a4)
        out.write(new_line)
        dih += 1

print('%d dihedrals' %(dih-1))