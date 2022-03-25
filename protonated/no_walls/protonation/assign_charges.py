
from pysimm import system, forcefield
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='lmps file to replace charges')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()

######################################################################################
# Get the charges using Gasteiger partial charges with GAFF in pysimm
s = system.read_lammps(args.lmps)
# s = system.read_lammps('protonated.lmps')
f = forcefield.Gaff()
s.forcefield = f.name

s.apply_charges(f, charges='gasteiger')

s.pair_style = 'lj/cut/coul/long 15.0'
s.bond_style = 'hybrid harmonic'
s.angle_style = 'hybrid harmonic'
s.dihedral_style = 'hybrid charmm multi/harmonic'

charges = {}
for i,p in enumerate(s.particles):
    charges[str(i+1)] = p.charge

######################################################################################
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

        if a_type == '14':
            charge = -1
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
            break

# Write all other sections unchanged
for line in f:
    out.write(line)

#############################################################################################
####################################### QUICK CHECKS ########################################
#############################################################################################

print('Total charge in the system is %.4f' %(total_charge))