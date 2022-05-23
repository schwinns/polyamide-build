import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def get_bonding(topology):

    top = open(topology, 'r')

    # Get to the atoms directive
    for line in top:
        if line.startswith('[ atoms ]'):
            break

    bonding = {}
    for line in top:
        if not line.startswith(';') and len(line.split()) > 0:
            atom = int(line.split()[0])
            bonding[atom] = {'bonded' : []}
        elif len(line.split()) == 0:
            break

    # Get to the bonding section
    for line in top:
        if line.startswith('[ bonds ]'):
            break

    for line in top:
        if not line.startswith(';') and len(line.split()) > 0:

            a1 = int(line.split()[0])
            a2 = int(line.split()[1])

            bonding[a1]['bonded'].append(a2)
            bonding[a2]['bonded'].append(a1)

        elif len(line.split()) == 0:
            break

    return bonding

################################# INPUTS ########################################

bulk_lims = np.array([2.0,6.0])     # bulk cutoffs in nm (z-direction)

frame_start = 391
frame_end = 401

traj = '../hydrate_center.xtc'             # input trajectory
gro = '../hydrate_center.gro'              # input coordinate file
gro_out = './gro_files/hydrate_center_'                   # output GRO filename

# Load trajectory
t = md.load(traj, top=gro)
top = t.topology
bulk_z = bulk_lims[1] - bulk_lims[0]

for f in np.arange(frame_start, frame_end):

    # Get atom indices of the bulk
    lb = np.where(t.xyz[f,:,2] > bulk_lims[0])[0]
    ub = np.where(t.xyz[f,:,2] < bulk_lims[1])[0]
    atom_idx = [i for i in lb if i in ub]

    # Fix any broken waters in the bulk
    for i in atom_idx:
        if top.atom(i).residue.is_water:
            res = top.atom(i).residue
            for a in res.atoms:
                if a.index not in atom_idx:
                    atom_idx.remove(i)
                    break

    # Define new unitcell vectors from bulk cutoffs
    unitcell_vectors = np.array([[t.unitcell_vectors[f,0,:],
                                  t.unitcell_vectors[f,1,:],
                                  [0,0,bulk_z]]])

    # Write GRO file for each frame
    t2 = t.atom_slice(atom_idx)[f]
    coords = t2.xyz
    coords[0,:,2] = coords[0,:,2] - bulk_lims[0]
    print('Writing gro file for frame %d out of %d to %s...' %(f, frame_end, gro_out + str(f) + '.gro') )
    with md.formats.GroTrajectoryFile(gro_out + str(f) + '.gro', 'w') as gro:
        gro.write(coords, t2.top, time=t[f].time, unitcell_vectors=unitcell_vectors)

