import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
from matplotlib import cm
from matplotlib import colors as mcolors
from matplotlib.colors import ListedColormap

myfolder = sys.argv[1]
framerate=10

cluster_data = h5py.File(myfolder + '/clusters.h5', 'r')
traj_data = h5py.File(myfolder + '/traj.h5', 'r')

times = np.array(cluster_data['/data/time'])
cluster_ids = np.array(cluster_data['/data/cluster_ids'])
pos = np.array(traj_data['/particles/all/position/value'])
nframes = times.shape[0]

cluster_data.close()
traj_data.close()

for t in range(nframes):

    if t%framerate==0:

        print(t)

        xcurr = pos[t,:,0]
        ycurr = pos[t,:,1]
        cids_curr = cluster_ids[t,:]

        #find clusters containing only one particle
        #unique, counts = np.unique(cids_curr, return_counts=True)
        #singletons = unique[counts==1]

        #cids_curr = np.where(np.in1d(cids_curr,singletons), 0, cids_curr)
        
        #define colormap
        maxnum = np.max(cids_curr)
        print(maxnum)
        mycolors = cm.get_cmap('Purples_r', maxnum)
        newcolors=mycolors(np.linspace(0,1,maxnum))
        #make the first few clusters distinct colors
        newcolors[0,:-1] = mcolors.to_rgb('Red')
        newcolors[1,:-1] = mcolors.to_rgb('Orange')
        newcolors[2,:-1] = mcolors.to_rgb('Yellow')
        newcolors[3,:-1] = mcolors.to_rgb('Green')
        newcolors[4,:-1] = mcolors.to_rgb('Blue')
        newcmp = ListedColormap(newcolors)

        #plot figure
        fig = plt.figure()
        plt.scatter(xcurr,ycurr,c = cids_curr,cmap=newcmp,s=0.2)
        plt.gca().set_aspect('equal')
        plt.xlim([-50,50])
        plt.ylim([-50,50])
        plt.savefig('frame_%04d.png' % (t//framerate),bbox_inches='tight',dpi=300)
        plt.close()