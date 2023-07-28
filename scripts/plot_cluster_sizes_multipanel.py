import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib import colors as mcolors
from matplotlib.colors import ListedColormap

potential='wca'
kT=0.0
phi=0.40
va=1.00
Lx=100.0
nx=200

rc1 = 1.010216
rc2 = 1.122462

taus=[0.01,0.10,0.30,1.00,10.00]
lambdas=[0.50,1.00,2.00,3.00,5.00,10.00]

nchunks = 10
colors = cm.viridis(np.linspace(0,1,nchunks))

fig, ax = plt.subplots(len(taus),len(lambdas),figsize=(10.0,7.5),sharex=True,sharey=True)

for i in range(len(taus)):
    for j in range(len(lambdas)):
        tau = taus[i]
        Lambda = lambdas[j]
        for k in range(1,nchunks):
            data = np.loadtxt('/scratch0/laynefrechette/active-noise-assembly-results/%s/kT=%.01f/phi=%.02f/va=%.02f/tau=%.02f/lambda=%.02f/Lx=%.01f_Ly=%.01f/nx=%d_ny=%d/seed=1/prod/cluster_hist_rc=%f_chunk=%d.txt' % (potential, kT, phi, va, tau, Lambda, Lx, Lx, nx, nx, rc2, k))
            if k==1:
                avg_data = data[:,1]
            else:
                avg_data += data[:,1]
        avg_data /=(nchunks-1)
        n_avg = np.dot(data[:,0],avg_data)
        n_avg_sq = np.dot(data[:,0]**2,avg_data)
        ax[i,j].plot(data[:,0],avg_data,color='blue')
        #ax[i,j].text(3*10**2, 0.1, r'$\langle n \rangle=%.1f$' % n_avg, ha='center')
        ax[i,j].text(2*10**2, 0.1, r'$\frac{\langle n^2 \rangle}{\langle n \rangle}=%.1f$' % (n_avg_sq/n_avg), ha='center',fontsize=10)

        ax[i,j].set_yscale('log')
        ax[i,j].set_xscale('log')

        if i==len(taus)-1:
            ax[i,j].set_xlabel(r'$n$')
        if j==0:
            ax[i,j].set_ylabel(r'$P(n)$')
        if i==0:
            ax[i,j].set_title(r'$\lambda=%.01f$' % Lambda)
        if j==len(lambdas)-1:
            ax2 = ax[i,j].twinx()
            ax2.set_ylabel(r'$\tau=%.02f$' % tau)
            ax2.set_yticks([])
        ax[i,j].set_xlim([1,5000])
#plt.legend()
plt.suptitle(r'CSDs for %s particles in active bath' % potential.upper())
#plt.tight_layout()
plt.savefig('cluster_hist_vary_lambda_tau.png',dpi=300,bbox_inches='tight')


#Compare different distance cutoffs
fig, ax = plt.subplots(len(taus),len(lambdas),figsize=(10.0,7.5),sharex=True,sharey=True)

for i in range(len(taus)):
    for j in range(len(lambdas)):
        tau = taus[i]
        Lambda = lambdas[j]
        for k in range(1,nchunks):
            data_rc1 = np.loadtxt('/scratch0/laynefrechette/active-noise-assembly-results/%s/kT=%.01f/phi=%.02f/va=%.02f/tau=%.02f/lambda=%.02f/Lx=%.01f_Ly=%.01f/nx=%d_ny=%d/seed=1/prod/cluster_hist_rc=%f_chunk=%d.txt' % (potential, kT, phi, va, tau, Lambda, Lx, Lx, nx, nx, rc1, k))
            if k==1:
                avg_data_rc1 = data_rc1[:,1]
            else:
                avg_data_rc1 += data_rc1[:,1]

            data_rc2 = np.loadtxt('/scratch0/laynefrechette/active-noise-assembly-results/%s/kT=%.01f/phi=%.02f/va=%.02f/tau=%.02f/lambda=%.02f/Lx=%.01f_Ly=%.01f/nx=%d_ny=%d/seed=1/prod/cluster_hist_rc=%f_chunk=%d.txt' % (potential, kT, phi, va, tau, Lambda, Lx, Lx, nx, nx, rc2, k))
            if k==1:
                avg_data_rc2 = data_rc2[:,1]
            else:
                avg_data_rc2 += data_rc2[:,1]

            data_rc3 = np.loadtxt('/scratch0/laynefrechette/active-noise-assembly-results/%s/kT=%.01f/phi=%.02f/va=%.02f/tau=%.02f/lambda=%.02f/Lx=%.01f_Ly=%.01f/nx=%d_ny=%d/seed=1/prod/cluster_hist_chunk=%d.txt' % (potential, kT, phi, va, tau, Lambda, Lx, Lx, nx, nx, k))
            if k==1:
                avg_data_rc3 = data_rc3[:,1]
            else:
                avg_data_rc3 += data_rc3[:,1]

        ax[i,j].plot(data[:,0],avg_data_rc1,color='red',label=r'$r_c=%f$' % rc1,zorder=10)
        ax[i,j].plot(data[:,0],avg_data_rc2,color='purple',label=r'$r_c=%f$' % rc2,zorder=9)
        ax[i,j].plot(data[:,0],avg_data_rc3,color='blue',label=r'$r_c=%f$' % (rc2*1.1),zorder=8)

        ax[i,j].set_yscale('log')
        ax[i,j].set_xscale('log')

        if i==len(taus)-1:
            ax[i,j].set_xlabel(r'$n$')
        if j==0:
            ax[i,j].set_ylabel(r'$P(n)$')
        if i==0:
            ax[i,j].set_title(r'$\lambda=%.01f$' % Lambda)
        if j==len(lambdas)-1:
            ax2 = ax[i,j].twinx()
            ax2.set_ylabel(r'$\tau=%.02f$' % tau)
            ax2.set_yticks([])
        ax[i,j].set_xlim([1,5000])
#plt.legend()
plt.suptitle(r'CSDs for %s particles in active bath' % potential.upper())
#plt.tight_layout()
plt.savefig('cluster_hist_vary_lambda_tau_compare_rc.png',dpi=300,bbox_inches='tight')