
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to replace charges')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

def check_N(N, atoms_top):

    # input the atom id for an N to check the other N
    #   if the other N is type NH, then MPD-T
    #   if the other N is type LN, then MPD-L

    # Positions on MPD ring
    #       _1_
    #     4    2
    #     |    |
    # N - 5_  _3 - N
    #       6

    n_NH = 0
    if atoms_top[N]['type'] == '8':
        n_NH += 1

    # check the other N on MPD to determine MPD-L/MPD-T    
    for a0 in atoms_top[N]['bonded']: 
        if atoms_top[a0]['type'] == '1': # 3/5 position
            for a1 in atoms_top[a0]['bonded']:
                if atoms_top[a1]['type'] == '1': # 2/4 and 6 positions
                    for a2 in atoms_top[a1]['bonded']:
                        if atoms_top[a2]['type'] == '1' and not a2 == a0: # 1 and 3/5 positions
                            for a3 in atoms_top[a2]['bonded']:
                                if atoms[a3]['type'] == '8':
                                    n_NH += 1
    
    return n_NH

def check_C(C, atoms_top):

    # input the atom id for a C to check the other Cs
    #   if n_CT = 0, then TMC-C
    #   if n_CT = 1, then TMC-L
    #   if n_CT = 2, then TMC-T

    # Positions on TMC ring
    #        C(14)
    #        |
    #       _1_
    #     4    2
    #     |    |
    # C - 5_  _3 - C
    #       6

    n_CT = 0
    if atoms_top[C]['type'] == '12':
        n_CT += 1
    
    for a in atoms_top[C]['bonded']: 
        if atoms_top[a]['type'] == '1': # find CA bonded to C (pos. 14) --> CA is position 1
            for a0 in atoms_top[a]['bonded']:
                if atoms_top[a0]['type'] == '1': # check CAs in 2/4
                    for a1 in atoms_top[a0]['bonded']:
                        if atoms_top[a1]['type'] == '1' and not a1 == a: # check CA in 3/5 but not 1
                            for a2 in atoms_top[a1]['bonded']:
                                if atoms_top[a2]['type'] == '12': # find carboxyl bonded to CA in 3/5
                                    n_CT += 1

    return n_CT

def find_N(atom, atoms_top):

    # input an atom id and find nearest N
    Ns = ['8','11']

    if atoms_top[atom]['type'] in Ns:
        return atom
    else:
        for a in atoms_top[atom]['bonded']: # check atoms bonded to atom
            if atoms_top[a]['type'] in Ns:
                return a

        for a in atoms_top[atom]['bonded']: # check atoms bonded to a level
            for a0 in atoms_top[a]['bonded']:
                if atoms_top[a0]['type'] in Ns:
                    return a0

        for a in atoms_top[atom]['bonded']: # check atoms bonded to a0 level
            for a0 in atoms_top[a]['bonded']:
                for a1 in atoms_top[a0]['bonded']:
                    if atoms_top[a1]['type'] in Ns:
                        return a1


def find_C(atom, atoms_top):

    # input an atom id and find nearest C
    Cs = ['10','12']

    if atoms_top[atom]['type'] in Cs:
        return atom
    else:
        for a in atoms_top[atom]['bonded']: # check atoms bonded to atom
            if atoms_top[a]['type'] in Cs:
                return a

        for a in atoms_top[atom]['bonded']: # check atoms bonded to a level
            for a0 in atoms_top[a]['bonded']:
                if atoms_top[a0]['type'] in Cs:
                    return a0

        for a in atoms_top[atom]['bonded']: # check atoms bonded to a0 level
            for a0 in atoms_top[a]['bonded']:
                for a1 in atoms_top[a0]['bonded']:
                    if atoms_top[a1]['type'] in Cs:
                        return a1

#######################################################################################
######################### CREATE DICTIONARIES WITH FULL TOPOLOGY ######################
#######################################################################################

f = open(args.lmps, 'r')
# f = open('2MPD-2TMC.lmps', 'r')

# Skip everything before atoms
for line in f:

    if line.startswith('Atoms'): # skip to atoms section
        break

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
            'bonds' : {}
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


f.close()

#######################################################################################
############################### ASSIGN PROPER CHARGES #################################
#######################################################################################

charges = {}
for atom in atoms:

    if atoms[atom]['type'] == '15':
        charges[atom] = -0.830

    elif atoms[atom]['type'] == '16':
        charges[atom] = 0.415

    elif atoms[atom]['type'] == '11': # LN

        n_NH = check_N(atom, atoms_top=atoms)
        if n_NH == 1:
            charges[atom] = -0.421
        elif n_NH == 0:
            charges[atom] = -0.4025
    
    elif atoms[atom]['type'] == '6': # HN

        for a in atoms[atom]['bonded']: 
            if atoms[a]['type'] == '8': # HN on NH
                charges[atom] = 0.3728
            elif atoms[a]['type'] == '11': # HN on LN in MPD-L

                # determine MPD-L/MPD-T
                n_NH = check_N(a, atoms_top=atoms)
                if n_NH == 1: # MPD-T
                    charges[atom] = 0.3577
                elif n_NH == 0: # MPD-L
                    charges[atom] = 0.3503

    elif atoms[atom]['type'] == '8': # NH
        
        charges[atom] = -0.9254

    elif atoms[atom]['type'] == '12': # CT

        n_CT = check_C(atom, atoms_top=atoms)
        if n_CT == 1:
            charges[atom] = 0.7387
        elif n_CT == 2:
            charges[atom] = 0.7369

    elif atoms[atom]['type'] == '10': # LC

        n_CT = check_C(atom, atoms_top=atoms)
        if n_CT == 0:
            charges[atom] = 0.50285
        elif n_CT == 1:
            charges[atom] = 0.5566
        elif n_CT == 2:
            charges[atom] = 0.5083

    elif atoms[atom]['type'] == '4': # O

        C = find_C(atom, atoms_top=atoms)
        n_CT = check_C(C, atoms_top=atoms)
        if n_CT == 0:
            charges[atom] = -0.46753
        elif n_CT == 1:
            
            if atoms[C]['type'] == '12': # if the C bonded to O is a CT
                charges[atom] = -0.5853
            else:
                charges[atom] = -0.5126

        elif n_CT == 2:

            if atoms[C]['type'] == '12': # if the C bonded to O is a CT
                charges[atom] = -0.5758
            else:
                charges[atom] = -0.4748

    ####### THIS ONE WILL CHANGE WHEN ACCOUNTING FOR DEPROTONATION ######
    elif atoms[atom]['type'] == '7': # OH

        C = find_C(atom, atoms_top=atoms)
        n_CT = check_C(C, atoms_top=atoms)
        if n_CT == 1:
            charges[atom] = -0.6148
        elif n_CT == 2:
            charges[atom] = -0.6373
    #####################################################################

    elif atoms[atom]['type'] == '9': # HO

        C = find_C(atom, atoms_top=atoms)
        n_CT = check_C(C, atoms_top=atoms)
        if n_CT == 1:
            charges[atom] = 0.4463
        elif n_CT == 2:
            charges[atom] = 0.4551

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
            N = bonded_types['11'][0]
            n_NH = check_N(N, atoms_top=atoms)
            if n_NH == 1:
                charges[atom] = 0.3883
            elif n_NH == 0:
                charges[atom] = -0.2175

        elif '10' in bonded_types: # CA bonded to LC
            C = find_C(atom, atoms_top=atoms)
            n_CT = check_C(C, atoms_top=atoms)
            if n_CT == 0:
                charges[atom] = -0.10058
            elif n_CT == 1:
                charges[atom] = -0.082
            elif n_CT == 2:
                charges[atom] = -0.0843

        elif '12' in bonded_types: # CA bonded to CT
            C = find_C(atom, atoms_top=atoms)
            n_CT = check_C(C, atoms_top=atoms)
            if n_CT == 1:
                charges[atom] = -0.0713
            elif n_CT == 2:
                charges[atom] = -0.0399

        elif '2' in bonded_types: # CA bonded to HA (CA, CA are other 2 bonds)
           
            bonded_types_next = {} # bonded_types_next[0 or 1] = [types of bonded atoms]
            i = 0
            for a in atoms[atom]['bonded']:
                if atoms[a]['type'] == '1': # check the types of the CAs on either side
                    bonded_types_next[i] = []
                    for a0 in atoms[a]['bonded']:
                        bonded_types_next[i].append(atoms[a0]['type'])
                    i += 1

            # Check the types of the atoms bonded to atom to determine charges
            #### MPD CAs ####
            if '2' in bonded_types_next[0] and '2' in bonded_types_next[1]: # MPD CA in 1 position
                N = find_N(atom, atoms_top=atoms)
                N_type = atoms[N]['type']
                if N_type == '8': # MPD-T
                    charges[atom] = 0.1293
                elif N_type == '11':
                    n_NH = check_N(N, atoms_top=atoms)
                    if n_NH == 0: # MPD-L
                        charges[atom] = -0.2418
                    elif n_NH == 1: # MPD-T
                        charges[atom] = 0.1293

            elif '11' in bonded_types_next[0] and '11' in bonded_types_next[1]: # MPD-L CA in 6 position
                charges[atom] = 0.1082

            elif '8' in bonded_types_next[0] and '11' in bonded_types_next[1]: # MPD-T CA in 6 position
                charges[atom] = -0.3154

            elif '8' in bonded_types_next[1] and '11' in bonded_types_next[0]: # MPD-T CA in 6 position
                charges[atom] = -0.3154

            elif '2' in bonded_types_next[0] and '11' in bonded_types_next[1]: # MPD CA in 2/4 position
                N = find_N(atom, atoms_top=atoms)
                N_type = atoms[N]['type']
                if N_type == '8': # MPD-T
                    charges[atom] = -0.2315
                elif N_type == '11':
                    n_NH = check_N(N, atoms_top=atoms)
                    if n_NH == 0: # MPD-L
                        charges[atom] = 0.3365
                    elif n_NH == 1: # MPD-T
                        charges[atom] = -0.2315

            elif '2' in bonded_types_next[1] and '11' in bonded_types_next[0]: # MPD CA in 2/4 position
                N = find_N(atom, atoms_top=atoms)
                N_type = atoms[N]['type']
                if N_type == '8': # MPD-T
                    charges[atom] = -0.2315
                elif N_type == '11':
                    n_NH = check_N(N, atoms_top=atoms)
                    if n_NH == 0: # MPD-L
                        charges[atom] = 0.3365
                    elif n_NH == 1: # MPD-T
                        charges[atom] = -0.2315

            elif '2' in bonded_types_next[0] and '8' in bonded_types_next[1]: # MPD-T CA in 2 position by NH
                charges[atom] = -0.2772

            elif '2' in bonded_types_next[1] and '8' in bonded_types_next[0]: # MPD-T CA in 2 position by NH
                charges[atom] = -0.2772

            #### TMC CAs ####
            elif '10' in bonded_types_next[0] and '10' in bonded_types_next[1]: # TMC CA between two CAs bonded to LCs
                C = find_C(atom, atoms_top=atoms)
                n_CT = check_C(C, atoms_top=atoms)
                if n_CT == 0:
                    charges[atom] = 0.06526
                elif n_CT == 1:
                    charges[atom] = 0.0423

            elif '12' in bonded_types_next[0] and '12' in bonded_types_next[1]: # TMC-T CA between two CAs bonded to CTs
                charges[atom] = 0.0344

            elif '10' in bonded_types_next[0] and '12' in bonded_types_next[1]: # TMC CA between CA bonded to LC and CA bonded to CT
                C = find_C(atom, atoms_top=atoms)
                n_CT = check_C(C, atoms_top=atoms)
                if n_CT == 1:
                    charges[atom] = 0.06005
                elif n_CT == 2:
                    charges[atom] = 0.0692

            elif '10' in bonded_types_next[1] and '12' in bonded_types_next[0]: # TMC CA between CA bonded to LC and CA bonded to CT
                C = find_C(atom, atoms_top=atoms)
                n_CT = check_C(C, atoms_top=atoms)
                if n_CT == 1:
                    charges[atom] = 0.06005
                elif n_CT == 2:
                    charges[atom] = 0.0692

                
    # print(charges)

#######################################################################################
###################### WRITE A NEW LAMMPS FILE WITH CORRECT CHARGES ###################
#######################################################################################

total_charge = 0
f = open(args.lmps, 'r')
out = open(args.output, 'w')
# Write Header and Coeffs
for line in f:

    if line.startswith('Atoms'):
        out.write(line)
        # out.write('\n')
        break

    else:
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

        if a_type == '2':
            charge = 0.0
        else:
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
            out.write(line)
            break

# Write all other sections unchanged
for line in f:
    out.write(line)

#############################################################################################
####################################### QUICK CHECKS ########################################
#############################################################################################

print('Total charge in the system is %.4f' %(total_charge))