
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lmps',
                    help='spce example file to relabel for PA system')
parser.add_argument('-o','--output',default='output.lmps',
                    help='output lmps filename')
args = parser.parse_args()


#######################################################################################
############################## WRITE A NEW LAMMPS FILE ################################
#######################################################################################

f = open(args.lmps, 'r')
out = open(args.output, 'w')
# Write Header and Coeffs
for line in f:

    if line.startswith('Atoms'):
        out.write(line)
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
        
        if a_type == '1':
            new_type = '15'
            charge = -0.83
        elif a_type == '2':
            new_type = '16'
            charge = 0.415
        
        x = l[4]
        y = l[5]
        z = l[6]

        new_line = " %s %s %s %.4f %s %s %s\n" %(a_id,mol_id,new_type,charge,x,y,z)
        out.write(new_line)

    else: # write blank lines
        out.write(line)

# Skip velocities if present
if vel:
    for line in f:
        if line.startswith('Bonds'):
            out.write(line)
            break


# Bonds section
for line in f:

    if line.startswith('Angles'):
        out.write(line)
        break

    elif len(line.split()) > 0:

        l = line.split()
        b_id = l[0]
        a1 = l[2]
        a2 = l[3]

        b_type = '19'

        new_line = " %s %s %s %s\n" %(b_id, b_type, a1, a2)
        out.write(new_line)

    else:
        out.write(line)

# Angles Section
for line in f:
    
    if len(line.split()) > 0:

        l = line.split()
        ang_id = l[0]
        a1 = l[2]
        a2 = l[3]
        a3 = l[4]

        ang_type = '24'

        new_line = " %s %s %s %s %s\n" %(ang_id, ang_type, a1, a2, a3)
        out.write(new_line)

    else:
        out.write(line)