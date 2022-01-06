
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to replace charges')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

def check_N(N, atoms_top=atoms):

    # input the atom id for an N to check the other N
    #   if the other N is type NH, then MPD-T
    #   if the other N is type LN, then MPD-L

    # Positions on MPD ring
    #       _1_
    #     4    2
    #     |    |
    # N - 5_  _3 - N
    #       6

    # check the other N on MPD to determine MPD-L/MPD-T
    for a0 in atoms_top[N]['bonded']: 
        if atoms_top[a0]['type'] == '1': # 3/5 position
            for a1 in atoms_top[a0]['bonded']:
                if atoms_top[a1]['type'] == '1': # 2/4 and 6 positions
                    for a2 in atoms_top[a1]['bonded']:
                        if atoms_top[a2]['type'] == '1': # 1 and 3/5 positions
                            for a3 in atoms_top[a2]['bonded']:
                                if atoms[a3]['type'] == '8': # MPD-T
                                    return 'MPD-T'
                                elif atoms[a3]['type'] == '11': # MPD-L
                                    return 'MPD-L'

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
            'monomer' : None
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
            'atoms' : [a1,a2]
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
            'atoms' : [a1,a2,a3]
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
            'atoms' : [a1,a2,a3,a4]
        }

#######################################################################################
############################### ASSIGN PROPER CHARGES #################################
#######################################################################################

# MPDs = {} # MPDs[MPD number][atom id] = charge
# m = 0
charges = {}
for atom in atoms:

    if atoms[atom]['type'] == '11': # LN

        charges[atom] = -0.4025
    
    elif atoms[atom]['type'] == '6': # HN

        for a in atoms[atom]['bonded']: 
            if atoms[a]['type'] == '8': # HN on NH
                charges[atom] = 0.3728 / 2 
            elif atoms[a]['type'] == '11': # HN on LN in MPD-L

                # determine MPD-L/MPD-T
                mono_type = check_N(a, atoms)
                if mono_type == 'MPD-T':
                    charges[atom] = 0.3577
                elif mono_type == 'MPD_L':
                    charges[atom] = 0.3503

    elif atoms[atom]['type'] == '1': # CA

        # save all the types bonded to each CA
        bonded_types = {}
        for a in atoms[atom]['bonded']:
            if atoms[a]['type'] not in bonded_types:
                bonded_types[atoms[a]['type']] = [a] # bonded_types[type] = [atom_id]
            else:
                bonded_types[atoms[a]['type']].append(a)

        # assign charges based on the bonded types
        if '8' in bonded_types: # CA bonded to NH
            charges[atom] = 0.5496

        elif '11' in bonded_types: # CA bonded to LN
            mono_type = check_N(bonded_types['11'][0],atoms)
            if mono_type == 'MPD-T':
                charges[atom] = 0.3883
            elif mono_type == 'MPD-L':
                charges[atom] = -0.2175

        elif 

        # m += 1
        # MPDs[m] = {atom : -0.4025} # add LN atom to MPD monomer

        # CAs = []
        # HAs = {} # HAs[CA] = HA bonded to that CA
        # for a0 in atoms[atom]['bonded']: # add CA and HN bonded to LN
        #     if atoms[a0]['type'] == '1':
        #         CAs.append(a0)
        #         MPDs[m] = {a0 : -0.2175}
        #     elif atoms[a0]['type'] == '6':
        #         HN1 = a0
        #         MPDs[m] = {a0 : 0.3503}

        # for a1 in atoms[CAs[0]]['bonded']: 
        #     CAs.append(a1)

        # for a2 in atoms[CAs[1]]['bonded']: # one of two CAs bonded to CAs[0]
        #     if atoms[a2]['type'] == '2':
        #         HAs[CAs[1]] = a2
        #     elif atoms[a2]['type'] == '1':
        #         CAs.append(a2)

        # for a2 in atoms[CAs[2]]['bonded']: # the other CA bonded to CAs[0]
        #     if atoms[a2]['type'] == '2':
        #         HAs[CAs[2]] = a2
        #     elif atoms[a2]['type'] == '1':
        #         CAs.append(a2)

        # for a2 in atoms[CAs[3]]['bonded']: # other CA bonded to CAs[1]
        #     if atoms[a2]['type'] == '2':
        #         HAs[CAs[1]] = a2
        #     elif atoms[a2]['type'] == '1':
        #         CAs.append(a2)

#######################################################################################
###################### WRITE A NEW LAMMPS FILE WITH CORRECT CHARGES ###################
#######################################################################################
# Write a new file with corrected charges
f = open(args.lmps,'r')
out = open(args.output, 'w')

total_charge = 0
# Write Header and Coeffs
for line in f:

    if line.startswith('Atoms'): # go to Atoms section
        out.write(line)
        break

    else: # everything before atoms
        out.write(line)

# Atoms section
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        vel = False
        out.write(line)
        break

    elif line.startswith('Velocities'):
        vel = True
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        a_id = l[0]
        mol_id = l[1]
        a_type = l[2]

        charge = charges[a_id]
        total_charge += charge
        
        x = l[4]
        y = l[5]
        z = l[6]

        new_line = " %s %s %s %.4f %s %s %s\n" %(a_id,mol_id,a_type,charge,x,y,z)
        out.write(new_line)

    else: # write blank lines
        out.write(line)

# Skip velocities if present
if vel:
    for line in f:
        if line.startswith('Bonds'):
            break

# Write all other sections unchanged
for line in f:
    out.write(line)

#############################################################################################
####################################### QUICK CHECKS ########################################
#############################################################################################

print('Total charge in the system is %.4f' %(total_charge))