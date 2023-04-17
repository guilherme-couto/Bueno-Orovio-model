import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import imageio.v2
import math

if len(sys.argv) != 2:
    print('Usage: python3 plot_ceq.py <dt>')
    sys.exit(1)

dt = sys.argv[1]

# Read data from files
t = []
timesfile = f'./simulation-files/sim-times-cable-eq-{dt}.txt'
with open(timesfile, 'r') as f:
    for line in f:
        t.append(float(line))

totalframes = len(t)

filename = f'./simulation-files/mm-cable-eq-{dt}.txt'
U = np.zeros((totalframes, 100))
with open(filename, 'r') as f:
    for n in range(totalframes):
        for i in range(100):
            U[n][i] = float(f.readline())


# Make plots
framerate = math.ceil(totalframes / 100)
plots = []
for n in range(totalframes):
    if n % framerate == 0:
        plotname = 'plot-' + str(n) + '.png'
        plots.append(plotname)
        
        plt.plot(U[n])
        plt.title(f'Cable Equation MM dt = {dt} t = {t[n]:.2f}')
        plt.xlabel('x')
        plt.ylabel('V (mV)')
        plt.ylim(-85, 50)
        
        plt.savefig(plotname)
        plt.close()

# Build gif
with imageio.v2.get_writer(f'./gif/gif-cable-eq-{dt}.gif', mode='I') as writer:
    for plot in plots:
        image = imageio.v2.imread(plot)
        writer.append_data(image)
        
# Remove files
for png in set(plots):
    os.remove(png)

    
    