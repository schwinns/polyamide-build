#!/usr/bin/env python

# Check total charge in lmps file

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to replace charges')                  
args = parser.parse_args()

#############################################################################################
####################################### USER INPUTS #########################################
#############################################################################################

# Keep track of the total charge in the system and how many of each atom type
total_charge = 0
n_types = {
    '1'  : 0, # CA
    '2'  : 0, # HA
    '3'  : 0, # C
    '4'  : 0, # O
    '5'  : 0, # N
    '6'  : 0, # HN
    '7'  : 0, # OH
    '8'  : 0, # NH
    '9'  : 0, # HO
    '10' : 0, # LC
    '11' : 0, # LN
    '12' : 0, # CT
}

#############################################################################################
################################### REWRITE WITH CHARGES ####################################
#############################################################################################

f = open(args.lmps,'r')

# Write Header and Coeffs
for line in f:

    if line.startswith('Atoms'): # go to Atoms section
        break

# Atoms section
for line in f:

    if line.startswith('Bonds'): # go to Bonds section
        vel = False
        break

    elif line.startswith('Velocities'):
        vel = True
        break

    elif len(line.split()) > 0: # if line is not blank

        l = line.split()

        a_type = l[2]
        charge = float(l[3])

        total_charge += charge
        n_types[a_type] += 1

#############################################################################################
####################################### QUICK CHECKS ########################################
#############################################################################################

print('Total charge in the system is %.4f' %(total_charge))
print('The number of each atom type:')
print(n_types)