import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
from matplotlib import cm
from matplotlib import colors as mcolors
from matplotlib.colors import ListedColormap

myfolder = sys.argv[1]

nchunks = 10
colors = cm.viridis(np.linspace(0,1,nchunks))

fig = plt.figure()
for i in range(nchunks):
    data = np.loadtxt(myfolder + '/cluster_hist_chunk=%d.txt' % i)
    plt.plot(data[:,0],data[:,1],color=colors[i],label='%d'%i)
plt.xlabel(r'$n$')
plt.ylabel(r'$P(n)$')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.savefig('cluster_hist.png',dpi=300)