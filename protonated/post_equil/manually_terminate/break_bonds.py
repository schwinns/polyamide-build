# Fix lmps file to break bonds across periodic boundaries

import argparse
import numpy as np
import random

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to move atoms')
parser.add_argument('-v','--verbose',action='store_true',
                    help='print completion percentages')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

f = open(args.lmps, 'r')
out = open(args.output, 'w')

#######################################################################################
############################ INITIAL INFORMATION FROM HEADER ##########################
#######################################################################################

# Get box lengths and write everything before atoms
print('Before breaking bonds:')
for line in f:

    if line.startswith('Atoms'):
        out.write(line)
        out.write('\n')
        break

    elif len(line.split()) == 2:
        out.write(line)
        new_line = line.split('\n')
        print('\t' + new_line[0])

    elif len(line.split()) == 4:

        l = line.split()
        lo = float(l[0])
        hi = float(l[1])

        box = hi - lo
        out.write(line)

    else:
        out.write(line)

print('\nBox dimensions (Angstroms):')
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

#######################################################################################
############################ LOCATE NEAREST LC-LN CROSSLINK BONDS #####################
#######################################################################################

for b in bonds:

    if args.verbose and int(b) % 100 == 0 or int(b) == len(bonds):
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
                if args.verbose:
                    print('Force break at iteration 50. Retrying at bond %s' %(b))
                break
            else:
                i += 1

        uncrosslink = b0
    

# Angle and dihedral types to change
to_change = {
    'angles' : {
        '15' : '22', # CA-LC-O  --> CA-CT-O
        '16' : '23', # CA-CA-LC --> CA-CA-CT
        '17' : '11', # LN-CA-CA --> NH-CA-CA
        '18' : '12' # CA-LN-HN --> CA-NH-HN
    },
    'dihedrals' : {
        '10' : '40', # CA-CA-CA-LC --> CA-CA-CA-CT
        '13' : '41', # LC-CA-CA-HA --> CT-CA-CA-HA
        '23' : '42', # CA-CA-CA-LC (imp) --> CA-CA-CA-CT (imp)
        '29' : '39', # CA-CA-LC-O --> CA-CA-CT-O
        '15' : '18', # LN-CA-CA-CA --> NH-CA-CA-CA
        '17' : '20', # LN-CA-CA-HA --> NH-CA-CA-HA
        '26' : '27', # CA-CA-CA-LN (imp) --> CA-CA-CA-NH (imp)
        '34' : '8' # HN-LN-CA-CA --> HN-NH-CA-CA
    }
}

# Mark the angles and dihedrals that should be deleted/changed in final lmps file
for bond in bonds:

    if bonds[bond]['delete']:

        a1 = bonds[bond]['atoms'][0]
        a2 = bonds[bond]['atoms'][1]

        # change LC and LN types to their terminated versions (CT, NH)
        if atoms[a1]['type'] == '10':
            CT = a1
            NH = a2
            atoms[CT]['type'] = '12'
            atoms[NH]['type'] = '8'
        elif atoms[a1]['type'] == '11':
            NH = a1
            CT = a2
            atoms[NH]['type'] = '8'
            atoms[CT]['type'] = '12'
       

        for angle in angles:
            if CT in angles[angle]['atoms'] and NH in angles[angle]['atoms']: # both in the angle, delete
                angles[angle]['delete'] = True
            elif CT in angles[angle]['atoms'] or NH in angles[angle]['atoms']: # if one in the angle, change
                old_type = angles[angle]['type']
                new_type = to_change['angles'][old_type]
                angles[angle]['type'] = new_type


        for dih in dihedrals:
            if CT in dihedrals[dih]['atoms'] and NH in dihedrals[dih]['atoms']:
                dihedrals[dih]['delete'] = True
            elif CT in dihedrals[dih]['atoms'] and NH in dihedrals[dih]['atoms']:
                old_type = dihedrals[dih]['type']
                new_type = to_change['dihedrals'][old_type]
                dihedrals[dih]['type'] = new_type

#######################################################################################
########### CREATE NEW ATOMS/BONDS/ANGLES/DIHEDRALS FOR TERMINATION GROUPS ############
#######################################################################################

# Get number of not broken bonds, angles, dihedrals
n_atoms = len(atoms)
n_bonds = 0
n_angles = 0
n_dihedrals = 0
n_broken = 0
for bond in bonds: 
    if not bonds[bond]['delete']:
        n_bonds += 1
    else:
        n_broken += 1

for angle in angles: 
    if not angles[angle]['delete']:
        n_angles += 1

for dihedral in dihedrals:
    if not dihedrals[dihedral]['delete']:
        n_dihedrals += 1

print('Broken crosslinks: %s' %(n_broken))

print('Before termination:')
print('\t%s atoms' %(n_atoms))
print('\t%s bonds' %(n_bonds))
print('\t%s angles' %(n_angles))
print('\t%s dihedrals\n' %(n_dihedrals))

# Place HN where broken LC is and OH where broken LN is
new_bonds = {}
for bond in bonds:

    if bonds[bond]['delete']:

        CA_TMC = [] # [CA-bonded-to-CT, CA-to-side, CA-to-side]
        CA_MPD = [] # [CA-bonded-to-NH, CA-to-side, CA-to-side]

        a1 = bonds[bond]['atoms'][0]
        a2 = bonds[bond]['atoms'][1]

        if atoms[a1]['type'] == '12':
            LC = a1
            LN = a2
        else:
            LC = a2
            LN = a1

        for a in atoms[LC]['bonded']: # get bonded atoms on TMC side for new angles, dihedrals
            if atoms[a]['type'] == '1':
                CA_TMC.append(a) # CA bonded to CT
                for a0 in atoms[a]['bonded']:
                    if atoms[a0]['type'] == '1':
                        CA_TMC.append(a0) # CAs on either side of CA-CT        
            elif atoms[a]['type'] == '4':
                O_TMC = a # CT=O

        for a in atoms[LN]['bonded']: # get bonded atoms on MPD side for new angles, dihedrals
            if atoms[a]['type'] == '1':
                CA_MPD.append(a) # CA bonded to NH
                for a0 in atoms[a]['bonded']:
                    if atoms[a0]['type'] == '1':
                        CA_MPD.append(a0) # CAs on either side of CA-NH     
            elif atoms[a]['type'] == '6':
                HN_MPD = a # HN-NH (original)

        # Generate neighbor lists for new HN and OH
        HN = dict(atoms[LC])
        OH = dict(atoms[LN])

        HN_neighbors = []
        for a in atoms:
            sq_dist = np.dot(atoms[a]['coords'] - HN['coords'],atoms[a]['coords'] - HN['coords'])
            if np.sqrt(sq_dist) < 5: # less than 5 angstroms away, it is a neighbor
                # print('Atom %s of type %s distance: %s' %(a, atoms[a]['type'], np.sqrt(sq_dist)))
                HN_neighbors.append(a)

        OH_neighbors = []
        for a in atoms:
            sq_dist = np.dot(atoms[a]['coords'] - OH['coords'],atoms[a]['coords'] - OH['coords'])
            if np.sqrt(sq_dist) < 5: # less than 5 angstroms away, it is a neighbor
                OH_neighbors.append(a)

        ####### HN ADDITION #######
        # Add HN to atoms
        HN['type'] = '6' # change type
        HN['bonded'] = [LN] # only bonded to broken LN
        HN['bonds'] = {str(n_bonds + 1) : LN}
    
        # Place the new atoms without overlapping
        close = True
        i = 0
        xyz = np.array([0,0,0])
        while close:

            i += 1
            min_dist = 100

            # check the distance from neighbors
            for a in HN_neighbors: 
                sq_dist = np.dot(atoms[a]['coords'] - HN['coords'],atoms[a]['coords'] - HN['coords'])
                
                if min_dist > np.sqrt(sq_dist): # save the smallest distance between neighbors
                    min_dist = np.sqrt(sq_dist)

                if min_dist < 2:
                    close = True
                else:
                    close = False

            # if i % 100 == 0:
            #             print('\t%d: %s --> %s Ang' %(i, HN['coords'], min_dist))

            if close: # if still close shift slightly in random direction
                xyz[0] = random.choice([-1,0,1])
                xyz[1] = random.choice([-1,0,1])
                xyz[2] = random.choice([-1,0,1])
                HN['coords'] = HN['coords'] + xyz*0.1
            
            elif args.verbose:
                print('Placed HN atom %d after %d iterations' %(n_atoms+1, i))
                break

            else:
                break

            if i == 5000:
                print('Could not place HN atom %d after %d iterations. Exiting...' %(n_atoms+1, i))
                exit()
            
        n_atoms += 1
        atoms[str(n_atoms)] = HN

        # Add HN-NH to bonds
        n_bonds += 1
        new_bonds[str(n_bonds)] = {
            'type' : '9',
            'atoms' : [LN,str(n_atoms)],
            'broken' : False,
            'delete' : False
        }

        # Add HN-NH-HN, CA-NH-HN to angles
        n_angles += 1
        angles[str(n_angles)] = { # HN-NH-HN
            'type' : '14',
            'atoms' : [LN,HN_MPD,str(n_atoms)],
            'delete' : False
        }
        n_angles += 1
        angles[str(n_angles)] = { # CA-NH-HN
            'type' : '12',
            'atoms' : [CA_MPD[0],LN,str(n_atoms)],
            'delete' : False
        }

        # Add CA-CA-NH-HN (x2), CA-NH-HN-HN (imp) to dihedrals
        n_dihedrals += 1
        dihedrals[str(n_dihedrals)] = { # CA-CA-NH-HN
            'type' : '8',
            'atoms' : [CA_MPD[1],CA_MPD[0],LN,str(n_atoms)],
            'delete' : False
        }
        n_dihedrals += 1
        dihedrals[str(n_dihedrals)] = { # CA-CA-NH-HN
            'type' : '8',
            'atoms' : [CA_MPD[2],CA_MPD[0],LN,str(n_atoms)],
            'delete' : False
        }
        n_dihedrals += 1
        dihedrals[str(n_dihedrals)] = { # CA-NH-HN-HN (imp)
            'type' : '28',
            'atoms' : [CA_MPD[0],LN,HN_MPD,str(n_atoms)],
            'delete' : False
        }

        ####### OH ADDITION #######
        OH['type'] = '7' # change type
        OH['bonded'] = [LC] # only bonded to broken LN
        OH['bonds'] = {str(n_bonds + 1) : LC}

        # Place the new atoms without overlapping
        close = True
        i = 0
        while close:

            i += 1
            min_dist = 100

            # check the distance from neighbors
            for a in OH_neighbors: 
                sq_dist = np.dot(atoms[a]['coords'] - OH['coords'],atoms[a]['coords'] - OH['coords'])
                
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
                OH['coords'] = OH['coords'] + xyz*0.1
            
            elif args.verbose:
                print('Placed OH atom %d after %d iterations' %(n_atoms+1, i))
                break

            else:
                break

            if i == 5000:
                print('Could not place OH atom %d after %d iterations. Exiting...' %(n_atoms+1, i))
                exit()

        n_atoms += 1
        atoms[str(n_atoms)] = OH

        # Add OH-CT to bonds
        n_bonds += 1
        new_bonds[str(n_bonds)] = {
            'type' : '14',
            'atoms' : [LC,str(n_atoms)],
            'broken' : False,
            'delete' : False
        }

        # Add O-CT-OH, CA-CT-OH to angles
        n_angles += 1
        angles[str(n_angles)] = { # O-CT-OH
            'type' : '21',
            'atoms' : [O_TMC,LC,str(n_atoms)],
            'delete' : False
        }
        n_angles += 1
        angles[str(n_angles)] = { # CA-CT-OH
            'type' : '19',
            'atoms' : [CA_TMC[0],LC,str(n_atoms)],
            'delete' : False
        }

        # Add CA-CA-CT-OH (x2), CA-CT-O-OH (imp) to dihedrals
        n_dihedrals += 1
        dihedrals[str(n_dihedrals)] = { # CA-CA-CT-OH
            'type' : '38',
            'atoms' : [CA_TMC[1],CA_TMC[0],LC,str(n_atoms)],
            'delete' : False
        }
        n_dihedrals += 1
        dihedrals[str(n_dihedrals)] = { # CA-CA-CT-OH
            'type' : '38',
            'atoms' : [CA_TMC[2],CA_TMC[0],LC,str(n_atoms)],
            'delete' : False
        }
        n_dihedrals += 1
        dihedrals[str(n_dihedrals)] = { # CA-CT-O-OH (imp)
            'type' : '21',
            'atoms' : [CA_MPD[0],LC,O_TMC,str(n_atoms)],
            'delete' : False
        }
   
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

print('After termination:')
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

        new_line = '%d %s %s %s\n' %(b, b_type, a1, a2)
        out.write(new_line)

for bond in new_bonds: # add the new bonds to the end

    b += 1
    b_type = new_bonds[bond]['type']
    a1 = new_bonds[bond]['atoms'][0]
    a2 = new_bonds[bond]['atoms'][1]

    new_line = '%d %s %s %s\n' %(b, b_type, a1, a2)
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

        new_line = '%d %s %s %s %s\n' %(ang, ang_type, a1, a2, a3)
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

        new_line = '%d %s %s %s %s %s\n' %(dih, dih_type, a1, a2, a3, a4)
        out.write(new_line)

print('\t%d dihedrals' %(dih))