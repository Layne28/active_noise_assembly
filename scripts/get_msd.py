#Measure MSD of particle in active noise flow

import numpy as np
import h5py
import sys
import numba

def main():
    myfolder = sys.argv[1]
    tmax = float(sys.argv[2])
    print('Computing msd for', myfolder)

    traj = h5py.File(myfolder + '/traj.h5')
    pos = np.array(traj['particles/all/position/value'])
    image = np.array(traj['particles/all/image/value'])
    edges = np.array(traj['particles/all/box/edges'])
    time = np.array(traj['particles/all/position/time'])
    traj.close()

    #unwrap the trajectory
    upos = unwrap(pos, image, edges)
    msd, data_times = get_msd(upos, time, tmax)

    print('Done.')

    np.savetxt(myfolder + '/msd.txt', np.c_[data_times, msd])

@numba.jit(nopython=True)
def get_msd(pos, time, tmax):

    dt = time[1]-time[0]
    frame_max = int(tmax/dt)
    data_times = np.arange(frame_max)*dt

    msd = np.zeros(frame_max)

    for t0 in range(time.shape[0]-frame_max):
        for t in range(frame_max):
            for d in range(2):
                #Careful here -- need to account for pbc
                msd[t] += (pos[t0+t][0][d]-pos[t0][0][d])*(pos[t0+t][0][d]-pos[t0][0][d])
    msd /= (time.shape[0]-frame_max)

    return msd, data_times

@numba.jit(nopython=True)
def unwrap(pos, image, edges):

    upos = np.zeros(pos.shape)
    for t in range(pos.shape[0]):
        for d in range(2):
            upos[t][0][d] = pos[t][0][d] + image[t][0][d]*edges[d]

    return upos

main()