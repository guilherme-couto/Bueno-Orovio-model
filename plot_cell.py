import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print('Usage: python3 plot_cell.py <dt>')
    sys.exit(1)

dt = sys.argv[1]

# Read data from files
t = []
timesfile = f'./simulation-files/sim-times-cell-adj-{dt}.txt'
with open(timesfile, 'r') as f:
    for line in f:
        t.append(float(line))

filename = f'./simulation-files/mm-cell-adj-{dt}.txt'
U = []
with open(filename, 'r') as f:
    for line in f:
        U.append(float(line))

# Make plot
plt.plot(t, U)
plt.title(f'Cell MM Adj')
plt.xlabel('t (ms)')
plt.ylabel('V (mV)')

plt.savefig('cell-adj.png')
plt.close()


    
    