#!/usr/bin/env python

# Fix Gromacs gro file to have proper atom numbers after packing with packmol

def write_gro_pos(pos, out, name='NA', box=None, ids=None, res=None, vel=None, ucell=None):
    """ write a .gro file from positions

    :param pos: xyz coordinates (natoms, 3)
    :param out: name of output .gro file
    :param name: name to give atoms being put in the .gro
    :param box: unitcell vectors. Length 9 list or length 3 list if box is cubic
    :param ids: name of each atom ordered by index (i.e. id 1 should correspond to atom 1)
    :param: res: name of residue for each atom
    :param: vel: velocity of each atom (natoms x 3 numpy array)
    :param: ucell: unit cell dimensions in mdtraj format (a 3x3 matrix)

    :type pos: np.ndarray
    :type out: str
    :type name: str
    :type box: list
    :type ids: list
    :type res: list
    :type vel: np.ndarray
    :type ucell: np.ndarray

    :return: A .gro file
    """

    if ucell is not None:
        box = [ucell[0, 0], ucell[1, 1], ucell[2, 2], ucell[0, 1], ucell[2, 0], ucell[1, 0], ucell[0, 2], ucell[1, 2],
               ucell[2, 0]]

    if box is None:  # to avoid mutable default
        box = [0., 0., 0.]

    with open(out, 'w') as f:

        f.write('This is a .gro file\n')
        f.write('%s\n' % pos.shape[0])

        for i in range(pos.shape[0]):
            if vel is not None:
                if ids is not None:
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n'.format((i + 1) % 100000, '%s' % name, '%s' % name,
                                                                            (i + 1) % 100000, pos[i, 0], pos[i, 1], pos[i, 2], vel[i, 0], vel[i, 1], vel[i, 2]))
                else:
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n'.format((i + 1) % 100000, '%s' % res[i], '%s' % ids[i],
                                                                            (i + 1) % 100000, pos[i, 0], pos[i, 1], pos[i, 2], vel[i, 0], vel[i, 1], vel[i, 2]))

            else:
                if ids is None:
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format((i + 1) % 100000, '%s' % name, '%s' % name,
                                                                            (i + 1) % 100000, pos[i, 0], pos[i, 1], pos[i, 2]))
                else:
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format((i + 1) % 100000, '%s' % res[i], '%s' % ids[i],
                                                                            (i + 1) % 100000, pos[i, 0], pos[i, 1], pos[i, 2]))
        for i in range(len(box)):
            f.write('{:10.5f}'.format(box[i]))

        f.write('\n')
        # f.write('{:10f}{:10f}{:10f}\n'.format(0, 0, 0))



import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-g','--gro',
                    help='gro file to rename atoms')
parser.add_argument('-o','--output',default='output.gro',
                    help='output gro filename')
args = parser.parse_args()

f = open(args.gro, 'r')
header = f.readline()
n_tot = int(f.readline().split()[0])

# out = open(args.output, 'w')
# out.write(header)
# out.write(n_tot)

xyz = np.zeros((n_tot,3))
box = [3.89842,3.89842,8.78931]
ids = []
res = []

# atoms = {} # atoms[atom_name] = # of atom type
number = 0
for line in f:
    
    if len(line.split()) != 3: # not box vectors

        l = line.split()
        molecule_name = l[0].strip('0123456789')
        # molecule_number = int(l[0].strip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
        atom_name = l[1].strip('0123456789')
        if len(l) == 6:
            x = float(l[3])
            y = float(l[4])
            z = float(l[5])
            xyz[number] = np.array([x,y,z])
        elif len(l) == 5:
            x = float(l[2])
            y = float(l[3])
            z = float(l[4])
            xyz[number] = np.array([x,y,z])
        else:
            print('Incorrect number of columns... gro may contain velocities')
            print('Line:%s' %(line) )
            exit()

        ids.append(atom_name)
        res.append(molecule_name)
        number += 1

        # if not atom_name in atoms:
        #     atoms[atom_name] = 0
        #     atom = atom_name
        # else:
        #     atoms[atom_name] += 1
        #     atom = atom_name + str(atoms[atom_name])

        # new_line = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(molecule_number % 100000, molecule_name, atom, number % 100000, x, y, z)
        # out.write(new_line)
        #out.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(molecule_number, molecule_name, atom, number, x, y, z) )

    # else:
        # out.write(line)

f.close()

write_gro_pos(xyz, args.output, box=box, ids=ids, res=res)