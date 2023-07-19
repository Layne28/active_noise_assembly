#Cluster displacement and velocity fields

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
    percentile_thresh = float(sys.argv[2])

    run_tests(in_folder)

    subfolder_list = []
    
    if(os.path.exists(in_folder)):
        for subfolder in sorted(sorted(glob.glob(in_folder + "/seed=*")),key=len):
            print(subfolder)
            subfolder_list.append(subfolder)
    else:
        print('Error: folder does not exist.')
        exit()

    #Load thresholds
    thresh_table = pd.read_csv(in_folder + '/thresh_disp_vel.txt', sep=' ')
    print(thresh_table)

    disp_thresh = thresh_table.loc[thresh_table['percentile']==percentile_thresh, 'disp'].item()
    vel_thresh = thresh_table.loc[thresh_table['percentile']==percentile_thresh, 'vel'].item()
    print(disp_thresh, vel_thresh)

    #Initialize data containers
    disp_cluster_sizes = [] #number of nodes in clusters
    disp_cluster_areas = [] #number of surface-exposed nodes in clusters
    all_num_disp_clusters = np.zeros(1)-1 #This will become bigger array. Set to -1 as check on whether it has been appended to yet
    vel_cluster_sizes = [] #number of nodes in clusters
    vel_cluster_areas = [] #number of surface-exposed nodes in clusters
    all_num_vel_clusters = np.zeros(1)-1 #This will become bigger array. Set to -1 as check on whether it has been appended to yet
    
    for subfolder in subfolder_list:
        print(subfolder)
        #Load data from reference lattice
        eq_traj = h5py.File(subfolder+'/equil/traj.h5', 'r')
        N, Lx, Ly, Lz, ref_lat = load_ref(eq_traj)
        bond_from_array = np.array(eq_traj['/parameters/vmd_structure/bond_from'])-1
        bond_to_array = np.array(eq_traj['/parameters/vmd_structure/bond_to'])-1
        adj_array = np.vstack((bond_from_array, bond_to_array)).T
        eq_traj.close()

        #Extract neighbor lists for each node
        #nbs = np.zeros((N,12), dtype=int) #assumes fixed connectivity of 12
        nbs = []
        for i in range(N):
            nbs_i1 = np.unique(adj_array[adj_array[:,0]==i][:,1])
            nbs_i2 = np.unique(adj_array[adj_array[:,1]==i][:,0])
            nbs_i = np.concatenate([nbs_i1,nbs_i2])
            #check that this is 12 for the present fcc lattice
            if nbs_i.shape[0]!=12:
                print('Warning: node does not have normal fcc connectivity.')
            nbs.append(nbs_i)
        if len(nbs)!=N:
            print('Warning: neighbor list length does not equal # of nodes.')

        #Load trajectory
        traj = h5py.File(subfolder+'/prod/traj.h5','r')

        #Get magnitude of displacements and velocities
        #use "unwrapped" positions to compute displacement field
        pos_unwrapped = get_pos_unwrapped(np.array(traj['/particles/all/position/value']), np.array(traj['/particles/all/image/value']), Lx, Ly, Lz)
        disp_all = pos_unwrapped-ref_lat
        disp_mag = la.norm(disp_all,axis=2)
        vel_all = np.array(traj['/particles/all/velocity/value'])
        vel_mag = la.norm(vel_all,axis=2)
        times = np.array(traj['/particles/all/velocity/time'])

        traj.close()

        #Threshold based on input percentile
        disp_threshed = np.where(disp_mag>disp_thresh,1,0)
        vel_threshed = np.where(vel_mag>vel_thresh,1,0)

        #Identify "clusters" as contiguous collections of high displacements/velocities
        traj_length = disp_threshed.shape[0]
        num_disp_clusters = np.zeros(traj_length)
        num_vel_clusters = np.zeros(traj_length)

        #Create files for dumping clusters
        uc_file = h5py.File(subfolder + '/prod/disp_clustered.h5', 'w')
        vc_file = h5py.File(subfolder + '/prod/vel_clustered.h5', 'w')
        uc_file.create_dataset("ref_lat", data=ref_lat)
        vc_file.create_dataset("ref_lat", data=ref_lat)

        for t in range(traj_length):
            if t%100==0:
                print('frame ', t)

            disp_cluster_id = sort_clusters(get_clusters(disp_threshed[t,:], nbs))
            num_disp_clusters[t] = np.max(disp_cluster_id)
            disp_unique, disp_counts = np.unique(disp_cluster_id, return_counts=True)

            vel_cluster_id = sort_clusters(get_clusters(vel_threshed[t,:], nbs))
            num_vel_clusters[t] = np.max(vel_cluster_id)
            vel_unique, vel_counts = np.unique(vel_cluster_id, return_counts=True)                

            #Record volume (# of nodes in cluster)
            #and (TODO) surface area (# of exposed nodes in cluster)
            disp_cluster_sizes.append(disp_counts[1:]) #exclude the "0" cluster
            vel_cluster_sizes.append(vel_counts[1:]) #exclude the "0" cluster

            #Save cluster-labeled data to file
            if t==0:
                uc_file.create_dataset('/data/time', data=np.array([times[t]]), chunks=True, maxshape=(None,))
                uc_file.create_dataset('/data/cluster_ids', data=np.array([disp_cluster_id]), chunks=True, maxshape=(None,N))
                uc_file.create_dataset('/data/value', data=np.array([disp_all[t,:,:]]), chunks=True, maxshape=(None,N,3))
                vc_file.create_dataset('/data/time', data=np.array([times[t]]), chunks=True, maxshape=(None,))
                vc_file.create_dataset('/data/cluster_ids', data=np.array([vel_cluster_id]), chunks=True, maxshape=(None,N))
                vc_file.create_dataset('/data/value', data=np.array([vel_all[t,:,:]]), chunks=True, maxshape=(None,N,3))
            else:
                uc_file['/data/time'].resize(uc_file['/data/time'].shape[0] + 1, axis=0)
                uc_file['/data/time'][-1] = times[t]
                uc_file['/data/cluster_ids'].resize(uc_file['/data/cluster_ids'].shape[0] + 1, axis=0)
                uc_file['/data/cluster_ids'][-1] = disp_cluster_id
                uc_file['/data/value'].resize(uc_file['/data/value'].shape[0] + 1, axis=0)
                uc_file['/data/value'][-1] = disp_all[t,:,:]
                vc_file['/data/time'].resize(vc_file['/data/time'].shape[0] + 1, axis=0)
                vc_file['/data/time'][-1] = times[t]
                vc_file['/data/cluster_ids'].resize(vc_file['/data/cluster_ids'].shape[0] + 1, axis=0)
                vc_file['/data/cluster_ids'][-1] = vel_cluster_id
                vc_file['/data/value'].resize(vc_file['/data/value'].shape[0] + 1, axis=0)
                vc_file['/data/value'][-1] = vel_all[t,:,:]

        if np.any(all_num_disp_clusters == -1):
            all_num_disp_clusters = num_disp_clusters
        else:
            all_num_disp_clusters = np.vstack([all_num_disp_clusters, num_disp_clusters])

        if np.any(all_num_vel_clusters == -1):
            all_num_vel_clusters = num_vel_clusters
        else:
            all_num_vel_clusters = np.vstack([all_num_vel_clusters, num_vel_clusters])

        uc_file.close()
        vc_file.close()

    #Get histograms
    disp_cluster_sizes = np.concatenate(disp_cluster_sizes)
    all_num_disp_clusters = all_num_disp_clusters.flatten()

    disp_size_hist, disp_size_bin_edges = np.histogram(disp_cluster_sizes, np.arange(0,N+2,1)-0.5, density=True)
    disp_size_bins = (disp_size_bin_edges[:-1]+disp_size_bin_edges[1:])/2

    disp_num_hist, disp_num_bin_edges = np.histogram(all_num_disp_clusters, np.arange(0,N+2,1)-0.5, density=True)

    vel_cluster_sizes = np.concatenate(vel_cluster_sizes)
    all_num_vel_clusters = all_num_vel_clusters.flatten()

    vel_size_hist, vel_size_bin_edges = np.histogram(vel_cluster_sizes, np.arange(0,N+2,1)-0.5, density=True)
    vel_size_bins = (vel_size_bin_edges[:-1]+vel_size_bin_edges[1:])/2

    vel_num_hist, vel_num_bin_edges = np.histogram(all_num_vel_clusters, np.arange(0,N+2,1)-0.5, density=True)

    #Write data
    np.savetxt(in_folder + '/cluster_hist_disp_thresh=%f.txt' % percentile_thresh, np.c_[disp_size_bins,disp_size_hist,disp_num_hist], header='bin size num')
    np.savetxt(in_folder + '/cluster_hist_vel_thresh=%f.txt' % percentile_thresh, np.c_[vel_size_bins,vel_size_hist,vel_num_hist], header='bin size num') 

def load_ref(traj):

    N = traj['/particles/all/position/value'].shape[1]
    Lx = traj['/particles/all/box/edges'][0]
    Ly = traj['/particles/all/box/edges'][1]
    Lz = traj['/particles/all/box/edges'][2]
    ref_lat = np.array(traj['/particles/all/position/value'][0,:])

    return N, Lx, Ly, Lz, ref_lat

def extract_frame(traj, t, Lx, Ly, Lz):
    
    pos = traj['/particles/all/position/value'][t,:]
    vel = traj['/particles/all/velocity/value'][t,:]
    tstep = traj['/particles/all/step'][t]
    time = traj['/particles/all/time'][t]
    
    return tstep, time, pos, vel

@numba.jit(nopython=True)
def get_pos_unwrapped(pos, image, Lx, Ly, Lz):

    shift = np.zeros(image.shape)
    shift[:,:,0] = image[:,:,0]*Lx
    shift[:,:,1] = image[:,:,1]*Ly
    shift[:,:,2] = image[:,:,2]*Lz

    pos_unwrapped = pos + shift
    
    return pos_unwrapped

@numba.jit(nopython=True)
def get_clusters(field, nbs):

    #Returns a numpy array specifying the index of the cluster
    #to which each node belongs

    N = field.shape[0]

    #print('Getting clusters...')
    cluster_id = np.zeros((N,),dtype=numba.int32)

    clusternumber = 0
    for i in range(N):
        if field[i]==1 and cluster_id[i]==0:
            clusternumber += 1
            cluster_id[i] = clusternumber
            harvest_cluster(clusternumber, i, cluster_id, field, nbs)

    return cluster_id

@numba.jit(nopython=True)
def harvest_cluster(clusternumber, ipart, cluster_id, field, nbs):

    for i in range(nbs[ipart].shape[0]):
        j = nbs[ipart][i]
        if field[j]==1 and cluster_id[j]==0:
            cluster_id[j] = clusternumber
            harvest_cluster(clusternumber, j, cluster_id, field, nbs)

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
        
def run_tests(in_folder):

    #Load data from reference lattice
    eq_traj = h5py.File(in_folder+'seed=1/equil/traj.h5','r')
    N, Lx, Ly, Lz, ref_lat = load_ref(eq_traj)
    bond_from_array = np.array(eq_traj['/parameters/vmd_structure/bond_from'])-1
    bond_to_array = np.array(eq_traj['/parameters/vmd_structure/bond_to'])-1
    adj_array = np.vstack((bond_from_array, bond_to_array)).T
    eq_traj.close()

    #Extract neighbor lists for each node
    #nbs = np.zeros((N,12), dtype=int) #assumes fixed connectivity of 12
    nbs = []
    for i in range(N):
        nbs_i1 = np.unique(adj_array[adj_array[:,0]==i][:,1])
        nbs_i2 = np.unique(adj_array[adj_array[:,1]==i][:,0])
        nbs_i = np.concatenate([nbs_i1,nbs_i2])
        #check that this is 12 for the present fcc lattice
        if nbs_i.shape[0]!=12:
            print('Warning: node does not have normal fcc connectivity.')
        nbs.append(nbs_i)
    if len(nbs)!=N:
        print('Warning: neighbor list length does not equal # of nodes.')

    all_zeros = np.zeros(N)
    all_ones = np.ones(N)

    cluster_id_zeros = get_clusters(all_zeros, nbs)
    cluster_id_ones = get_clusters(all_ones, nbs)

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