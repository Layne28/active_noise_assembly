#Get Cluster Size Distributions

import numpy as np
import numpy.linalg as la
import pylab as plt
import sys
import os
import glob
import h5py
import pandas as pd
import numba
import math

def main():

    in_folder = sys.argv[1]
    nchunks = 10 #make nchunks histograms for sequential trajectory segments

    #Initialize data containers
    cluster_sizes = [] #number of nodes in clusters
    cluster_areas = [] #number of surface-exposed nodes in clusters
    all_num_clusters = np.zeros(1)-1 #This will become bigger array. Set to -1 as check on whether it has been appended to yet
    
    #Create file for dumping clusters
    #TODO: change this to modifying traj.h5
    cluster_file = h5py.File(in_folder + '/clusters.h5', 'r')
    times = np.array(cluster_file['/data/time'])
    cluster_ids = np.array(cluster_file['/data/cluster_ids'])
    cluster_file.close()

    traj_length = times.shape[0]
    N = cluster_ids.shape[1]
    num_clusters = np.zeros(traj_length)
    cluster_sizes = []

    for t in range(traj_length):
        if t%100==0:
            print('frame ', t)

        cluster_id = cluster_ids[t,:]
        num_clusters[t] = np.max(cluster_id)
        unique, counts = np.unique(cluster_id, return_counts=True) #counts gives cluster size distribution at time t
        #print(counts)

        #Record volume (# of particles in cluster)
        #and (TODO) surface area (# of exposed particles in cluster)
        cluster_sizes.append(counts)

    #Get histograms
    hist_list = []
    size_list = []
    chunk_size = int(traj_length/nchunks)
    for n in range(nchunks):
        
        cluster_size_hist, size_bin_edges = np.histogram(cluster_sizes[n], np.arange(0,N+2,1)-0.5, density=True)
        print(cluster_size_hist)
        size_bins = (size_bin_edges[:-1]+size_bin_edges[1:])/2
        print(size_bins)

        #Write data
        np.savetxt(in_folder + '/cluster_hist_chunk=%d.txt' % n, np.c_[size_bins,cluster_size_hist], header='bin size')

main()