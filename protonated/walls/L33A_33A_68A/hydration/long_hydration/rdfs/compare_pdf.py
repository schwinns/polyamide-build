# Script to quickly plot .xvg files from GROMACS
# Uses GromacsWrapper so should be in LLC-env or gro_wrap

import numpy as np
from gromacs.formats import XVG
import mdtraj as md
import json
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from time import time

# compute_rdf taken from MDTraj and modified to use truncated unit cell
def compute_rdf(traj, pairs, r_range=None, bin_width=0.005, n_bins=None,
                periodic=True, opt=True, bulk_lims=None, scale=True, scaling_factors=None):
    """Compute radial distribution functions for pairs in every frame.
    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute radial distribution function in.
    pairs : array-like, shape=(n_pairs, 2), dtype=int
        Each row gives the indices of two atoms.
    r_range : array-like, shape=(2,), optional, default=(0.0, 1.0)
        Minimum and maximum radii.
    bin_width : float, optional, default=0.005
        Width of the bins in nanometers.
    n_bins : int, optional, default=None
        The number of bins. If specified, this will override the `bin_width`
         parameter.
    periodic : bool, default=True
        If `periodic` is True and the trajectory contains unitcell
        information, we will compute distances under the minimum image
        convention.
    opt : bool, default=True
        Use an optimized native library to compute the pair wise distances.
    Returns
    -------
    r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radii values corresponding to the centers of the bins.
    g_r : np.ndarray, shape=(np.diff(r_range) / bin_width - 1), dtype=float
        Radial distribution function values at r.
    See also
    --------
    Topology.select_pairs
    """
    if r_range is None:
        r_range = np.array([0.0, 1.0])
    r_range = md.utils.ensure_type(r_range, dtype=np.float64, ndim=1, name='r_range',
                          shape=(2,), warn_on_cast=False)
    if n_bins is not None:
        n_bins = int(n_bins)
        if n_bins <= 0:
            raise ValueError('`n_bins` must be a positive integer')
    else:
        n_bins = int((r_range[1] - r_range[0]) / bin_width)

    distances = md.geometry.distance.compute_distances(traj, pairs, periodic=periodic, opt=opt)
    print('\tFinished computing distances... distances array.shape = (%d, %d)' %(distances.shape[0], distances.shape[1]))
    if scale:
        g_r, edges = np.histogram(distances, range=r_range, bins=n_bins, weights=scaling_factors)
        r = 0.5 * (edges[1:] + edges[:-1])
        print('\tFinished the weighted histogram')
    else:
        g_r, edges = np.histogram(distances, range=r_range, bins=n_bins)
        r = 0.5 * (edges[1:] + edges[:-1])

    if bulk_lims is not None:
        unitcell_vol = traj.unitcell_volumes * ((bulk_lims[1] - bulk_lims[0]) / traj.unitcell_lengths[:,2])
    else:
        print('\tDefaulting to scaling by full simulation cell...')
        unitcell_vol = traj.unitcell_volumes

    # Normalize by volume of the spherical shell.
    # See discussion https://github.com/mdtraj/mdtraj/pull/724. There might be
    # a less biased way to accomplish this. The conclusion was that this could
    # be interesting to try, but is likely not hugely consequential. This method
    # of doing the calculations matches the implementation in other packages like
    # AmberTools' cpptraj and gromacs g_rdf.
    V = (4 / 3) * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
    norm = len(pairs) * np.sum(1.0 / unitcell_vol) * V
    g_r = g_r.astype(np.float64) / norm  # From int64.
    return r, g_r


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
            bonding[atom] = {'bonded' : [],
                             'one_away' : [],
                             'two_away' : [],
                             'exclusions' : []}
        elif len(line.split()) == 0:
            break

    # Get to the pairs section
    for line in top:
        if line.startswith('[ pairs ]'):
            break

    for line in top:
        if not line.startswith(';') and len(line.split()) > 0:

            a1 = int(line.split()[0])
            a2 = int(line.split()[1])

            bonding[a1]['two_away'].append(a2)
            bonding[a2]['two_away'].append(a1)
            bonding[a1]['exclusions'].append(a2)
            bonding[a2]['exclusions'].append(a1)

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
            bonding[a1]['exclusions'].append(a2)
            bonding[a2]['exclusions'].append(a1)

        elif len(line.split()) == 0:
            break

    # Add 1-3 pairs
    for atom in bonding:

        for a1 in bonding[atom]['bonded']:
            for a2 in bonding[atom]['bonded']:
                
                if not a1 == a2:
                    bonding[a1]['one_away'].append(a2)
                    bonding[a1]['exclusions'].append(a2)

    return bonding


#################################################################################
############################## SIMULATION RDFS ##################################
#################################################################################

################################# INPUTS ########################################

excl = True                         # if True, do not calculate interatomic distances for 1-2, 1-3, 1-4 atoms in molecules
water = True                        # if True, include water in RDF
bulk = True                         # if True, only calculate RDF for the bulk defined by bulk_lims (only consider atoms within cutoff in last frame)
bulk_lims = np.array([1.5,4.5])     # bulk cutoffs in nm (z-direction)
scale = True                        # if True, scale the RDF by atomic form factors

frame_start = 400                   # first frame to calculate RDF
frame_end = 401                     # last frame to calculate RDF
timing = True                       # if True, display timing information
plot = True                         # if True, show final RDF plot

traj = '../hydrate.xtc'             # input trajectory
gro = '../hydrate.gro'              # input coordinate file
topology = '../PA_hydrated.top'     # input PA topology
json_factors = './form_factors.json'# json file with form factors
filename = './output.xvg'           # output RDF filename

#################################################################################

# Load trajectory
start = time()
print('\n\n--------------------------- PROGRESS ---------------------------')
print("Loading trajectory '%s' with topology '%s'" %(traj, gro) )
t = md.load(traj, top=gro)
top = t.topology

# Get necessary information for the input parameters
if excl:
    print("Excluding bonded atoms. Getting bonding information from '%s'" %(topology))
    bonding = get_bonding(topology)
if bulk:
    lb = np.where(t.xyz[-1,:,2] > bulk_lims[0])[0]
    ub = np.where(t.xyz[-1,:,2] < bulk_lims[1])[0]

    atom_idx = np.array([i for i in lb if i in ub])
    print('Calculating RDF for bulk membrane. Found %d atoms within %.2f nm to %.2f nm' %(len(atom_idx), bulk_lims[0], bulk_lims[1]) )

else:
    atom_idx = [atom.index for atom in top.atoms]

# Calculate the average form factor
if scale:
    print('Finding pairs of atoms and pre-calculating form factors...')
    form_factors = json.load(open(json_factors))
    factors = np.zeros(len(atom_idx))
    for i, idx in enumerate(atom_idx):
        f_i = form_factors[top.atom(idx).element.symbol]
        factors[i] = f_i

    avg_f = factors.mean()

else:
     print('Finding pairs of atoms...')
     scaling_factors = None

load_time = time()
# Apply the inputs to select only desired atom pairs
scaling = []
pairs = []
for i in atom_idx:
    for j in atom_idx:
        if i != j:
            f_i = form_factors[top.atom(i).element.symbol]
            f_j = form_factors[top.atom(j).element.symbol]

            if excl:

                if top.atom(i).residue.is_water or top.atom(j).residue.is_water: # if either atom is water
                    if top.atom(i).residue.resSeq != top.atom(j).residue.resSeq: # check if they are same water molecule
                        pairs.append([i,j])
                        scaling.append(f_i*f_j)
                        
                else: # if not water, check if they should be excluded
                    if j not in bonding[i]['exclusions']:
                        pairs.append([i,j])
                        scaling.append(f_i*f_j)
            else:
                pairs.append([i,j])
                scaling.append(f_i*f_j)

scaling_factors = np.zeros((frame_end - frame_start, len(pairs)))
for f in range(frame_end - frame_start):
    scaling_factors[f,:] = np.array(scaling) / avg_f**2

pair_time = time()
print('Computing the RDF from frame %d (%.1f ps) to frame %d (%.1f ps)...' %(frame_start, t[frame_start-1].time, frame_end, t[frame_end-1].time) )
r, g_r = compute_rdf(t[frame_start:frame_end], pairs, bulk_lims=bulk_lims, scale=scale, scaling_factors=scaling_factors)

print("Writing the RDF data to '%s'" %(filename))
data = np.vstack((r*10, g_r))
xvg = XVG(array=data, names=['r [Å]', 'g(r)'])
xvg.write(filename)
print('----------------------------------------------------------------')

rdf_time = time()
if timing:
    print('\n\n----------------------- TIMING BREAKDOWN -----------------------')
    print('\tLoading trajectory and inputs:\t\t%f' %(load_time - start))
    print('\tFinding pairs:\t\t\t\t%f' %(pair_time - load_time))
    print('\tRDF computations:\t\t\t%f' %(rdf_time - pair_time))
    print('\n\tTotal time:\t\t\t\t%f' %(rdf_time - start))
    print('----------------------------------------------------------------')

if plot:
    fig, ax = plt.subplots(1,1)
    plt.plot(r*10, g_r, c='b', label='Simulation')
    plt.xlim(0,10)
    # plt.ylim(-5,10)
    plt.ylabel('g(r)')
    plt.xlabel('r (Å)')

    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.25))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.axhline(c='k', lw=0.5)

#################################################################################
############################## EXPERIMENTAL RDFS ################################
#################################################################################

################################# INPUTS ########################################

filename = './TN4'   # input file with G(r) data
rho = 0.0953         # atomic number density
rho_sim = True       # if True, calculate the number density from simulation

#################################################################################

# Calculate the number density from the atom indices in the bulk simulation in final frame
if rho_sim:
    xy = t.unitcell_lengths[-1,0]*10
    z = (bulk_lims[1] - bulk_lims[0])*10
    n_atoms = len(atom_idx)

    rho = n_atoms / xy**2 / z
    print('\nNumber density of bulk from simulation: %.4f atoms/Å^3' %(rho) )

# Read in data
exp_data = np.loadtxt(filename, comments='#')

# Clean data
min_idx = np.where(exp_data[:,0] == np.min(exp_data[:,0]))[0]
max_idx = np.where(exp_data[:,0] == np.max(exp_data[:,0]))[0]

## Data is from min_idx[0] to max_idx[0]
## Base line is from min_idx[1] to max_idx[1]

clean = np.zeros((max_idx[0] - min_idx[0], 2))
clean[:,0] = exp_data[min_idx[0]:max_idx[0],0]
clean[:,1] = exp_data[min_idx[0]:max_idx[0],1] / (4*np.pi*clean[:,0]*rho) + 1

# Plot the data
if plot:
    plt.plot(clean[:,0], clean[:,1], c='r', label='Experiment') # custom labels
    plt.legend()
    plt.show()