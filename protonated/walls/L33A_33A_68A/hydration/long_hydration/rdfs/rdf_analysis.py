# RDF analysis with MDTraj
# Scales by f_i(0)*f_j(0) / <f>^2

import matplotlib.pyplot as plt
import gromacs as gro
from matplotlib.ticker import MultipleLocator
import mdtraj as md
from mdtraj.utils import ensure_type
from mdtraj.geometry.distance import compute_distances
import numpy as np
import json

# compute_rdf taken from MDTraj and modified to scale distances
def compute_rdf(traj, pairs, r_range=None, bin_width=0.005, n_bins=None,
                periodic=True, opt=True):
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
    r_range = ensure_type(r_range, dtype=np.float64, ndim=1, name='r_range',
                          shape=(2,), warn_on_cast=False)
    if n_bins is not None:
        n_bins = int(n_bins)
        if n_bins <= 0:
            raise ValueError('`n_bins` must be a positive integer')
    else:
        n_bins = int((r_range[1] - r_range[0]) / bin_width)

    distances = compute_distances(traj, pairs, periodic=periodic, opt=opt)
    scaled_distances = scale_distances(traj, pairs, distances, json_factors='./form_factors.json')
    g_r, edges = np.histogram(scaled_distances, range=r_range, bins=n_bins)
    r = 0.5 * (edges[1:] + edges[:-1])

    # Normalize by volume of the spherical shell.
    # See discussion https://github.com/mdtraj/mdtraj/pull/724. There might be
    # a less biased way to accomplish this. The conclusion was that this could
    # be interesting to try, but is likely not hugely consequential. This method
    # of doing the calculations matches the implementation in other packages like
    # AmberTools' cpptraj and gromacs g_rdf.
    V = (4 / 3) * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
    norm = len(pairs) * np.sum(1.0 / traj.unitcell_volumes) * V
    g_r = g_r.astype(np.float64) / norm  # From int64.
    return r, g_r


def scale_distances(traj, pairs, distances, json_factors='./form_factors.json'):
    
    """Scale radial distribution distances for pairs in every frame.
    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute radial distribution function in.
    pairs : array-like, shape=(n_pairs, 2), dtype=int
        Each row gives the indices of two atoms.
    distances : np.ndarray, shape=(n_frames, n_pairs), dtype=float
        The distance, in each frame, between each pair of atoms
    json_factors : filename, dtype=string
        Name of the json file containing f(0) atomic form factor for all atoms, default=form_factors.json
    Returns
    -------
    scaled_distances : np.ndarray, shape=(n_frames, n_pairs), dtype=float
        The distance, in each frame, between each pair of atoms, scaled by atomic form factor
    """

    form_factors = json.load(open(json_factors))
    top = traj.topology
    
    scaled_distances = np.zeros(distances.shape)
    elements = {}
    indices = []
    for p, pair in enumerate(pairs):
        atom_i = top.atom(pair[0])
        atom_j = top.atom(pair[1])
        
        element_i = atom_i.element.symbol
        element_j = atom_j.element.symbol
        
        f_i = form_factors[element_i]
        f_j = form_factors[element_j]
        
        scaled_distances[:,p] = distances[:,p]*f_i*f_j
        
        if element_i not in elements and pair[0] not in indices:
            elements[element_i] = 1
            indices.append(pair[0])
        elif pair[0] not in indices:
            elements[element_i] += 1
            indices.append(pair[0])
        if element_j not in elements and pair[1] not in indices:
            elements[element_j] = 1
            indices.append(pair[1])
        elif pair[1] not in indices:
            elements[element_j] += 1
            indices.append(pair[1])

    avg_f = 0
    tot_elem = 0
    for element in elements:
        avg_f += elements[element]*form_factors[element]
        tot_elem += elements[element]
    
    avg_f = avg_f / tot_elem
    
    return scaled_distances / avg_f**2


def get_bonding(topology):

    top = open(topology, 'r')

    # Get to the atoms directive
    for line in top:
        if line.startswith('[ atoms ]'):
            break

    for line in top:
        if not line.startswith(';') and len(line.split()) > 0:
            atom = int(line.split()[0])
            bonding[atom] = {'bonded' : [],
                             'one_away': [],
                             'two_away' : [],
                             'three_away' : []}

        elif len(line.split()) == 0:
            break

    # Get to the bonding section
    for line in top:
        if line.startswith('[ bonds ]'):
            break

    bonding = {}
    for line in top:
        if not line.startswith(';') and len(line.split()) > 0:

            a1 = int(line.split()[0])
            a2 = int(line.split()[1])

            bonding[a1]['bonded'].append(a2)
            bonding[a2]['bonded'].append(a1)



#################################################################################
########################### MAIN USE OF FUNCTIONS ###############################
#################################################################################

t = md.load('../hydrate.xtc', top='../hydrate.gro')
top = t.topology

atom_idx = [atom.index for atom in top.atoms]

# Get bonding environment for each atom in PA membrane (saving 3 levels of bonding info)
bonding = get_bonding('../PA_hydrated.top')


# print(bonding)

# pairs = []

# for i in atom_idx:
#     for j in atom_idx:
#         if i != j:
#             pairs.append([i,j])

# r, g_r = compute_rdf(t, pairs)

# fig, ax = plt.subplots(1,3, figsize=(18,5))

# ax[0].set_title('OO')
# ax[0].plot(dataO[0,:]*10, dataO[1,:], label='Unscaled')
# ax[0].plot(r_OO*10, g_OO, label='Scaled')

# ax[1].set_title('OH')
# ax[1].plot(dataO[0,:]*10, dataO[4,:], label='Unscaled')
# ax[1].plot(r_OH*10, g_OH, label='Scaled')

# ax[2].set_title('HH')
# ax[2].plot(dataH[0,:]*10, dataH[4,:], label='Unscaled')
# ax[2].plot(r_HH*10, g_HH, label='Scaled')
                
# for i in range(3):
#     ax[i].set_xlim(0,10)
# #     ax[i].set_ylim(-0.1,3)
#     ax[i].set_xlabel('r [Å]')
#     ax[i].set_ylabel('g(r)')
#     ax[i].legend();