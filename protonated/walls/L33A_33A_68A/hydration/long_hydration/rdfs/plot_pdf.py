# Script to quickly plot .xvg files from GROMACS
# Uses GromacsWrapper so should be in LLC-env or gro_wrap

import numpy as np
import gromacs as gro
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-xvg','--xvg',
                help='rdf file to plot against experimental data')
parser.add_argument('-f','--file',
                    help='txt file with experimental data')
args = parser.parse_args()

# Import the .xvg file
xvg = gro.fileformats.XVG(args.xvg)
sim_data = xvg.array
r = sim_data[0,:] #*10
rdf = sim_data[1,:] #- sim_data[1,np.where(r == 10)[0]]

fig, ax = plt.subplots(1,1)
plt.plot(r, rdf, c='b', label='Simulation')
plt.xlim(0,10)
# plt.ylim(-0.25,5)
plt.ylabel('g(r)')
plt.xlabel('r (A)')

ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.25))
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(0.25))
ax.axhline(c='k', lw=0.5)

# plt.savefig(args.xvg.split('.')[0] + '.png')

# Read in data
exp_data = np.loadtxt(args.file, comments='#')

# Clean data
min_idx = np.where(exp_data[:,0] == np.min(exp_data[:,0]))[0]
max_idx = np.where(exp_data[:,0] == np.max(exp_data[:,0]))[0]

## Data is from min_idx[0] to max_idx[0]
## Base line is from min_idx[1] to max_idx[1]

# number density?
# rho = 0.0953
rho = 0.1087

# shift_idx = np.where(exp_data[:,0] == 1)[0]
# clean = np.zeros((max_idx[0] - shift_idx[0], 2))
# clean[:,0] = exp_data[shift_idx[0]:max_idx[0],0]
# clean[:,1] = exp_data[shift_idx[0]:max_idx[0],1] / (4*np.pi*clean[:,0]*rho) + 1

clean = np.zeros((max_idx[0] - min_idx[0], 2))
clean[:,0] = exp_data[min_idx[0]:max_idx[0],0]
clean[:,1] = exp_data[min_idx[0]:max_idx[0],1] / (4*np.pi*clean[:,0]*rho) + 1

# Plot the data
plt.plot(clean[:,0], clean[:,1], c='r', label='Experiment') # custom labels
plt.legend()

# plt.savefig(args.file.split('.')[0] + '_sim_compare.png')

plt.show()


