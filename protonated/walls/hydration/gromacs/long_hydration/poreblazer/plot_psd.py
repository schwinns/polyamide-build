#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-p','--psd',
                    help='txt file with pore size distribution from PoreBlazer')
parser.add_argument('-o','--output',default='output.png',
                    help='output png filename')
args = parser.parse_args()

psd = open(args.psd,'r')

# Get all data
x = []
y = []
for line in psd:

    if not line.split()[0] == '#': # ignore comments

        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))

x = np.array(x) # convert to numpy arrays
y = np.array(y)

# Plot the data
fig = plt.figure()
plt.plot(x/10,y)
plt.xlabel('distance (nm)')
plt.ylabel('PSD')
plt.savefig(args.output)