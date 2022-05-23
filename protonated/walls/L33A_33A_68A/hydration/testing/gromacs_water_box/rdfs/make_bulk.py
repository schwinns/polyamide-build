import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt


################################# INPUTS ########################################

bulk_lims = np.array([1.0,2.0])     # bulk cutoffs in nm (z-direction)

frame_start = 0
frame_end = 51

traj = '../npt.xtc'             # input trajectory
gro = '../npt.gro'              # input coordinate file
gro_out = './gro_files/output_'                   # output GRO filename

######################### WRITE BULK GRO FILES ################################

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
    print('Writing gro file for frame %d out of %d to %s...' %(f, frame_end-1, gro_out + str(f) + '.gro') )
    with md.formats.GroTrajectoryFile(gro_out + str(f) + '.gro', 'w') as gro:
        gro.write(coords, t2.top, time=t[f].time, unitcell_vectors=unitcell_vectors)

