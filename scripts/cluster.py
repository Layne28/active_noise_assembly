#Cluster particles

import numpy as np
import numpy.linalg as la
import pylab as plt
import sys
import os
import glob
import h5py
import pandas as pd
import numba

def main():

    in_folder = sys.argv[1]
    #rc = float(sys.argv[2]) #default: 2^(1/6) plus a small number
    rc=pow(2,1.0/6)*1.2

    #run_tests(in_folder)

    subfolder_list = [in_folder]
    
    '''
    if(os.path.exists(in_folder)):
        for subfolder in sorted(sorted(glob.glob(in_folder + "/seed=*")),key=len):
            print(subfolder)
            subfolder_list.append(subfolder)
    else:
        print('Error: folder does not exist.')
        exit()
    '''

    #Initialize data containers
    cluster_sizes = [] #number of nodes in clusters
    cluster_areas = [] #number of surface-exposed nodes in clusters
    all_num_clusters = np.zeros(1)-1 #This will become bigger array. Set to -1 as check on whether it has been appended to yet
    
    for subfolder in subfolder_list:
        print(subfolder)

        #Load trajectory
        traj = h5py.File(subfolder+'/prod/traj.h5','r')

        N = traj['/particles/all/position/value'].shape[1]
        Lx = traj['/particles/all/box/edges'][0]
        Ly = traj['/particles/all/box/edges'][1]
        edges = np.array([Lx,Ly])

        pos= np.array(traj['/particles/all/position/value'])
        times = np.array(traj['/particles/all/position/time'])

        traj.close()

        traj_length = times.shape[0]
        num_clusters = np.zeros(traj_length)

        #Create file for dumping clusters
        cluster_file = h5py.File(subfolder + '/prod/clusters.h5', 'w')

        for t in range(traj_length):
            if t%100==0:
                print('frame ', t)

            cluster_id = sort_clusters(get_clusters(pos[t,:,:], edges, rc))
            num_clusters[t] = np.max(cluster_id)
            unique, counts = np.unique(cluster_id, return_counts=True)          

            #Record volume (# of particles in cluster)
            #and (TODO) surface area (# of exposed particles in cluster)
            cluster_sizes.append(counts[1:]) #exclude the "0" cluster

            #Save cluster-labeled data to file
            if t==0:
                cluster_file.create_dataset('/data/time', data=np.array([times[t]]), chunks=True, maxshape=(None,))
                cluster_file.create_dataset('/data/cluster_ids', data=np.array([cluster_id]), chunks=True, maxshape=(None,N))
                cluster_file.create_dataset('/data/value', data=np.array([pos[t,:,:]]), chunks=True, maxshape=(None,N,3))
            else:
                cluster_file['/data/time'].resize(cluster_file['/data/time'].shape[0] + 1, axis=0)
                cluster_file['/data/time'][-1] = times[t]
                cluster_file['/data/cluster_ids'].resize(cluster_file['/data/cluster_ids'].shape[0] + 1, axis=0)
                cluster_file['/data/cluster_ids'][-1] = cluster_id
                cluster_file['/data/value'].resize(cluster_file['/data/value'].shape[0] + 1, axis=0)
                cluster_file['/data/value'][-1] = pos[t,:,:]

        if np.any(all_num_clusters == -1):
            all_num_clusters = num_clusters
        else:
            all_num_clusters = np.vstack([all_num_clusters, num_clusters])

        cluster_file.close()

    #Get histograms
    cluster_sizes = np.concatenate(cluster_sizes)
    all_num_clusters = all_num_clusters.flatten()

    cluster_size_hist, size_bin_edges = np.histogram(cluster_sizes, np.arange(0,N+2,1)-0.5, density=True)
    size_bins = (size_bin_edges[:-1]+size_bin_edges[1:])/2

    num_hist, num_bin_edges = np.histogram(all_num_clusters, np.arange(0,N+2,1)-0.5, density=True)

    #Write data
    np.savetxt(in_folder + '/cluster_hist.txt', np.c_[size_bins,cluster_size_hist,num_hist], header='bin size num')

@numba.jit(nopython=True)
def get_clusters(pos, edges, rc):

    #Returns a numpy array specifying the index of the cluster
    #to which each node belongs

    N = pos.shape[0]

    #print('Getting clusters...')
    cluster_id = np.zeros((N,),dtype=numba.int32)

    clusternumber = 0
    for i in range(N):
        if cluster_id[i]==0:
            clusternumber += 1
            cluster_id[i] = clusternumber
            harvest_cluster(clusternumber, i, cluster_id, pos, edges, rc)

    return cluster_id

@numba.jit(nopython=True)
def harvest_cluster(clusternumber, ipart, cluster_id, pos, edges, rc):

    pos1 = np.array([pos[ipart,0],pos[ipart,1]]) #pos[ipart,:]
    for j in range(pos.shape[0]): #this is slow, could make faster with neighbor list
        pos2 = np.array([pos[j,0],pos[j,1]])#pos[j,:]
        #rij = get_min_dist(pos[ipart,:],pos[j,:], edges[0], edges[1])
        rij = get_min_dist(pos1,pos2, edges[0], edges[1])
        if j!=ipart and rij<=rc and cluster_id[j]==0:
            cluster_id[j] = clusternumber
            harvest_cluster(clusternumber, j, cluster_id, pos, edges, rc)

def sort_clusters(cluster_id_arr):

    #Re-order cluster ids from largest to smallest
    sizes = np.zeros((np.max(cluster_id_arr)+1,2),dtype=int)
    for i in range(sizes.shape[0]):
        sizes[i,0] = i
        sizes[i,1] = cluster_id_arr[cluster_id_arr==i].shape[0]
    sorted_sizes = sizes[1:]
    sorted_sizes = sorted_sizes[(-sorted_sizes[:,1]).argsort()] #sort in descending order, excluding zero cluster
    size_map = {} #make map from cluster id to size rank
    #print('sizes:', sorted_sizes)
    for i in range(sorted_sizes.shape[0]):
        size_map[sorted_sizes[i,0]] = i+1
    size_map[0] = 0

    #Now rename cluster ids
    cluster_id_sorted = np.array([size_map[a] for a in cluster_id_arr])
    if sizes[1:,1].shape[0]>0 and sizes[1:,1].max() != cluster_id_sorted[cluster_id_sorted==1].shape[0]:
        print('Error in re-labeling clusters! Biggest cluster does not have cluster id = 1.')
        exit()

    return cluster_id_sorted
        
@numba.jit(nopython=True) 
def get_min_disp(r1, r2, Lx, Ly):
    arr1 = np.array([Lx/2,Ly/2])
    arr2 = np.array([-Lx/2,-Ly/2])
    arr3 = np.array([Lx, Ly])
    rdiff = r1-r2
    rdiff = np.where(rdiff>arr1, rdiff-arr3, rdiff)
    rdiff = np.where(rdiff<arr2, rdiff+arr3, rdiff)
    return rdiff

@numba.jit(nopython=True) 
def get_min_dist(r1, r2, Lx, Ly):
    rdiff = get_min_disp(r1,r2,Lx,Ly)
    return la.norm(rdiff)

def run_tests(in_folder):

    N = 100

    all_zeros = np.zeros(N)
    all_ones = np.ones(N)

    edges=np.array([10,10])
    rc=1

    cluster_id_zeros = get_clusters(all_zeros,edges,rc)
    cluster_id_ones = get_clusters(all_ones,edges,rc)

    if(np.max(cluster_id_zeros)!=0):
        print('Test failed! Configuration with no nodes above threshold has cluster.')
        exit()
    if(np.max(cluster_id_ones)!=1):
        print('Test failed! Configuration with all nodes above threshold has no or more than one cluster.')
        exit()
    
    unique, counts = np.unique(cluster_id_ones, return_counts=True)
    if(counts[0]!=N):
        print('Test failed! Configuration with all nodes does not have cluster with N nodes.')
        exit()

    print('Tests passed!')
    
main()