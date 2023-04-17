import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import imageio.v2
import math

if len(sys.argv) != 2:
    print('Usage: python3 plot_cell.py <dt>')
    sys.exit(1)

dt = sys.argv[1]

# Read data from files
t = []
timesfile = f'./simulation-files/sim-times-cell-{dt}.txt'
with open(timesfile, 'r') as f:
    for line in f:
        t.append(float(line))

filename = f'./simulation-files/mm-cell-{dt}.txt'
U = []
with open(filename, 'r') as f:
    for line in f:
        U.append(float(line))

# Make plot
plt.plot(t, U)
plt.title(f'Cell MM')
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')

plt.savefig('Cell.png')
plt.close()


    
    