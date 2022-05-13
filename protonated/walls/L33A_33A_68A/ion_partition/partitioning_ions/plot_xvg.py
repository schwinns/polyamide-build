# Script to quickly plot .xvg files from GROMACS
# Uses GromacsWrapper so should be in LLC-env or gro_wrap

# Use as:
#           python plot_xvg.py <filename> -x xmin xmax -y ymin ymax

import numpy as np
import gromacs as gro
import sys
import matplotlib.pyplot as plt

# Get the .xvg file to plot from the command line
filename = sys.argv[1]

# Import the .xvg file
xvg = gro.fileformats.XVG(filename)
# xmax = xvg.array[0][-1]
data = xvg.array

# Plot using the GromacsWrapper and save as 'filename'.png
fig = plt.figure()
# plot = xvg.plot(color='Dark2', maxpoints=None)
plt.plot(data[0,:],data[1,:] / 6693, label='PA')
plt.plot(data[0,:],data[2,:] / 5, label='Na')
plt.plot(data[0,:],data[3,:] / 5, label='Cl')
plt.plot(data[0,:],data[4,:] / 7542, label='Water')

# plot.set_xlim(xmax-50000, xmax)

# name = filename[0:-4]
# plot.set_ylabel(name)
# plot.set_xlabel('time (ps)')
# plt.savefig(name + '_last50ns.png')

plt.ylabel('Number density')
plt.xlabel('z coordinate (nm)')
plt.legend()
plt.show()
