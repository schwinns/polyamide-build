# Fix lmps file to break bonds across periodic boundaries

import argparse
import numpy as np
import random

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps input file')
parser.add_argument('-p','--proton',default=0.5,type=float,
                    help='percentage of carboxyl groups that need protonation')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

f = open(args.lmps, 'r')
out = open(args.output, 'w')

#######################################################################################
############################ INITIAL INFORMATION FROM HEADER ##########################
#######################################################################################

# Write everything before atoms
for line in f:

    if line.startswith('Atoms'):
        out.write(line)
        out.write('\n')
        break

    else:
        out.write(line)

#######################################################################################
######################### CREATE DICTIONARIES WITH FULL TOPOLOGY ######################
#######################################################################################

# Create dictionaries with full topology
atoms = {}
for line in f: # build dictionary of all atoms

    if line.startswith('Velocities'): # end of atoms section
        vel = True
        break

    elif line.startswith('Bonds'):
        vel = False
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
            'bonds' : {},
            'delete' : False
        }

# skip velocities
if vel:
    for line in f:
        if line.startswith('Bonds'):
            break

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
############################### PROTONATE CARBOXYL GROUPS #############################
#######################################################################################

# First, get carboxyls that are already protonated
protonated = []
unprotonated = []
for atom in atoms:

    if atoms[atom]['type'] == '7': # if it is an OH 

        unprotonated.append(atom)

        for a in atoms[atom]['bonded']: # check the atoms bonded to the OH

            if atoms[a]['type'] == '9' and atom not in protonated: # if there is a hydrogen bonded
                unprotonated.remove(atom)
                protonated.append(atom)
        

# print('%d oxygens are protonated' %(len(protonated)))
# print('%d oxygens are unprotonated' %(len(unprotonated)))
# print('%.3f of the oxygens are protonated' %(len(protonated) / (len(protonated) + len(unprotonated)) ))

n_atoms = len(atoms)
n_bonds = len(bonds)
n_angles = len(angles)
n_dihedrals = len(dihedrals)

# Protonate the unprotonated oxygens
for atom in unprotonated:

    protonated.append(atom)
    
    # get atoms bonded for new angles, dihedrals
    CT = atoms[atom]['bonded'][0]
    for a in atoms[CT]['bonded']:
        if atoms[a]['type'] == '1':
            CA = a
        elif atoms[a]['type'] == '4':
            O = a

    HO = dict(atoms[atom])
    HO_neighbors = []
    for a in atoms:
        sq_dist = np.dot(atoms[a]['coords'] - HO['coords'],atoms[a]['coords'] - HO['coords'])
        if np.sqrt(sq_dist) < 5: # less than 5 angstroms away, it is a neighbor
            HO_neighbors.append(a)

    HO['type'] = '9'
    HO['bonded'] = [atom]
    HO['bonds'] = {str(n_bonds+1) : atom}

    # Place the new atoms without overlapping
    close = True
    i = 0
    xyz = np.array([0,0,0])
    while close:

        i += 1
        min_dist = 100

        # check the distance from neighbors
        for a in HO_neighbors: 
            sq_dist = np.dot(atoms[a]['coords'] - HO['coords'],atoms[a]['coords'] - HO['coords'])
            
            if min_dist > np.sqrt(sq_dist): # save the smallest distance between neighbors
                min_dist = np.sqrt(sq_dist)

            if min_dist < 2:
                close = True
            else:
                close = False

        if close: # if still close shift slightly in random direction
            xyz[0] = random.choice([-1,0,1])
            xyz[1] = random.choice([-1,0,1])
            xyz[2] = random.choice([-1,0,1])
            HO['coords'] = HO['coords'] + xyz*0.1

        else:
            break

        if i == 5000:
            print('Could not place HO atom %d after %d iterations. Exiting...' %(n_atoms+1, i))
            exit()
            
    n_atoms += 1
    atoms[str(n_atoms)] = HO

    # Add new HO-OH bond
    n_bonds += 1
    bonds[str(n_bonds)] = {
        'type' : '15',
        'atoms' : [atom,str(n_atoms)],
        'broken' : False,
        'delete' : False
    }

    # Add CT-OH-HO to angles
    n_angles += 1
    angles[str(n_angles)] = {
        'type' : '20',
        'atoms' : [CT,atom,str(n_atoms)],
        'delete' : False
    }

    # Add CA-CT-OH-HO, O-CT-OH-HO to dihedrals
    n_dihedrals += 1
    dihedrals[str(n_dihedrals)] = { # CA-CT-OH-HO
        'type' : '35',
        'atoms' : [CA,CT,atom,str(n_atoms)],
        'delete' : False
    }
    n_dihedrals += 1
    dihedrals[str(n_dihedrals)] = { # O-CT-OH-HO
        'type' : '36',
        'atoms' : [O,CT,atom,str(n_atoms)],
        'delete' : False
    }

#######################################################################################
############################## DEPROTONATE CARBOXYL GROUPS ############################
#######################################################################################

print('Before deprotonation:')
print('\t%d oxygens are protonated' %(len(protonated)))
print('\t%d oxygens are unprotonated' %(len(unprotonated)))
print('\t%.3f of the oxygens are protonated\n' %(len(protonated) / (len(protonated) + len(unprotonated)) ))

percent_prot = 1.0
unprotonated = []
while abs(percent_prot - args.proton) > 0.01:

    p = random.randint(0, len(protonated) - 1) # select a random oxygen to deprotonate
    OH = protonated[p]
    protonated.remove(OH)
    unprotonated.append(OH)
    # atoms[OH]['type'] = '14' # change type for charge assignments in next step

    percent_prot = len(protonated) / (len(protonated) + len(unprotonated))

    for a in atoms[OH]['bonded']:
        if atoms[a]['type'] == '9':
            HO = a
            atoms[a]['delete'] = True

    for bond in bonds:
        if HO in bonds[bond]['atoms']:
            bonds[bond]['delete'] = True

    for angle in angles:
        if HO in angles[angle]['atoms']:
            angles[angle]['delete'] = True

    for dih in dihedrals:
        if HO in dihedrals[dih]['atoms']:
            dihedrals[dih]['delete'] = True


print('After deprotonation:')
print('\t%d oxygens are protonated' %(len(protonated)))
print('\t%d oxygens are unprotonated' %(len(unprotonated)))
print('\t%.3f of the oxygens are protonated\n' %(len(protonated) / (len(protonated) + len(unprotonated)) ))


#######################################################################################
#################### WRITE A NEW LAMMPS FILE WITH CORRECT PROTONATION #################
#######################################################################################

# Header and coefficient were written above

# Atoms section
a = 0
atom_ids = {} # atom_ids[old_id] = new_id
for atom in atoms:

    if not atoms[atom]['delete']:

        a += 1
        atom_ids[atom] = a
        a_id = a
        mol_id = atoms[atom]['molecule']
        a_type = atoms[atom]['type']
        charge = atoms[atom]['charge']
        x = atoms[atom]['coords'][0]
        y = atoms[atom]['coords'][1]
        z = atoms[atom]['coords'][2]

        new_line = '%s %s %s %s %f %f %f\n' %(a_id, mol_id, a_type, charge, x, y, z)
        out.write(new_line)

print('\t%s atoms' %(a))

# Bonds section
out.write('\nBonds\n\n')
b = 0
for bond in bonds: # need to renumber bonds because some are deleted

    if not bonds[bond]['delete']:
    
        b += 1
        b_type = bonds[bond]['type']
        a1 = bonds[bond]['atoms'][0]
        a2 = bonds[bond]['atoms'][1]

        new_a1 = atom_ids[a1]
        new_a2 = atom_ids[a2]

        new_line = '%d %s %s %s\n' %(b, b_type, new_a1, new_a2)
        out.write(new_line)


print('\t%s bonds' % (b))

# Angles section
out.write('\nAngles\n\n')
ang = 0
for angle in angles: # need to renumber angles because some are deleted

    if not angles[angle]['delete']:
    
        ang += 1
        ang_type = angles[angle]['type']
        a1 = angles[angle]['atoms'][0]
        a2 = angles[angle]['atoms'][1]
        a3 = angles[angle]['atoms'][2]

        new_a1 = atom_ids[a1]
        new_a2 = atom_ids[a2]
        new_a3 = atom_ids[a3]

        new_line = '%d %s %s %s %s\n' %(ang, ang_type, new_a1, new_a2, new_a3)
        out.write(new_line)

print('\t%d angles' %(ang))

# Dihedrals section
out.write('\nDihedrals\n\n')
dih = 0
for dihedral in dihedrals: # need to renumber dihedrals because some are deleted

    if not dihedrals[dihedral]['delete']:
    
        dih += 1
        dih_type = dihedrals[dihedral]['type']
        a1 = dihedrals[dihedral]['atoms'][0]
        a2 = dihedrals[dihedral]['atoms'][1]
        a3 = dihedrals[dihedral]['atoms'][2]
        a4 = dihedrals[dihedral]['atoms'][3]

        new_a1 = atom_ids[a1]
        new_a2 = atom_ids[a2]
        new_a3 = atom_ids[a3]
        new_a4 = atom_ids[a4]

        new_line = '%d %s %s %s %s %s\n' %(dih, dih_type, new_a1, new_a2, new_a3, new_a4)
        out.write(new_line)

print('\t%d dihedrals' %(dih))