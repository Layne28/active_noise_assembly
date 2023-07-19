#Create conf file given parameters

import argparse
import os

def main():

    parser = argparse.ArgumentParser(description='Write active noise conf file with given input parameters.')
    parser.add_argument('confdir',
                        help='Directory where conf file will be written.')
    #System args
    parser.add_argument('--rho', default='1.414213562373095')
    parser.add_argument('--a', default='1.414213562373095')
    parser.add_argument('--Nx', default='5')
    parser.add_argument('--Ny', default='5')
    parser.add_argument('--Nz', default='5')
    parser.add_argument('--Lax', default='5')
    parser.add_argument('--Lay', default='5')
    parser.add_argument('--Laz', default='5')
    parser.add_argument('--is_p_x', default='1')
    parser.add_argument('--is_p_y', default='1')
    parser.add_argument('--is_p_z', default='1')
    parser.add_argument('--node_protocol', default='fcc')
    parser.add_argument('--spring_protocol', default='uniform')
    parser.add_argument('--l0', default='1.0')
    parser.add_argument('--K', default='1.0')

    #Solver args
    parser.add_argument('--dt', default='0.002')
    parser.add_argument('--gamma', default='1.0')
    parser.add_argument('--va', default='1.0')

    #Observer args
    parser.add_argument('--output_dir', default='active-network-results')
    parser.add_argument('--network_freq', default='100')
    parser.add_argument('--thermo_freq', default='100')

    #LabBench args
    parser.add_argument('--experiment', default='standard')
    parser.add_argument('--equil_steps', default='50000')
    parser.add_argument('--production_steps', default='100000')
    parser.add_argument('--info_freq', default='10000')

    #Active Noise Generator args
    parser.add_argument('--nx', default='32')
    parser.add_argument('--ny', default='32')
    parser.add_argument('--nz', default='32')
    parser.add_argument('--tau', default='1.0')
    parser.add_argument('--Lambda', default='1.0')

    args = parser.parse_args()
    print(args.__dict__)

    try:
        os.makedirs(args.confdir)
        print('Made directory.')
    except OSError as e:
        print('Directory exists')

    print('Writing conf file...')
    with open(args.confdir + '/active_net_va=%f_tau=%f_lambda=%f_Nx=%d_Ny=%d_Nz=%d_nx=%d.conf' % (float(args.va),float(args.tau),float(args.Lambda), int(args.Nx), int(args.Ny), int(args.Nz), int(args.nx)), 'w') as f:

        f.write('#System\n')
        f.write('rho = %s\n' % args.rho)
        f.write('a = %s\n' % args.a)
        f.write('Nx = %s\n' % args.Nx)
        f.write('Ny = %s\n' % args.Ny)
        f.write('Nz = %s\n' % args.Nz)
        f.write('Lax = %s\n' % args.Lax)
        f.write('Lay = %s\n' % args.Lay)
        f.write('Laz = %s\n' % args.Laz)
        f.write('is_p_x = %s\n' % args.is_p_x)
        f.write('is_p_y = %s\n' % args.is_p_y)
        f.write('is_p_z = %s\n' % args.is_p_z)
        f.write('node_protocol = %s\n' % args.node_protocol)
        f.write('spring_protocol = %s\n' % args.spring_protocol)
        f.write('l0 = %s\n' % args.l0)
        f.write('K = %s\n' % args.K)
        f.write('\n')

        f.write('#Solver\n')
        f.write('dt = %s\n' % args.dt)
        f.write('gamma = %s\n' % args.gamma)
        f.write('va = %s\n' % args.va)
        f.write('\n')

        f.write('#Observer\n')
        f.write('output_dir = %s\n' % args.output_dir)
        f.write('network_freq = %s\n' % args.network_freq)
        f.write('thermo_freq = %s\n' % args.thermo_freq)
        f.write('\n')

        f.write('#LabBench\n')
        f.write('experiment = %s\n' % args.experiment)
        f.write('equil_steps = %s\n' % args.equil_steps)
        f.write('production_steps = %s\n' % args.production_steps)
        f.write('info_freq = %s\n' % args.info_freq)
        f.write('\n')

        f.write('#Active Noise Generator\n')
        f.write('nx = %s\n' % args.nx)
        f.write('ny = %s\n' % args.ny)
        f.write('nz = %s\n' % args.nz)
        f.write('tau = %s\n' % args.tau)
        f.write('lambda = %s\n' % args.Lambda)

if __name__=='__main__':
    main()
