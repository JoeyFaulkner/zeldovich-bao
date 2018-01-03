import zeldovich as Z
import cic_dens_wrapper
import numpy as N
import matplotlib.pyplot as plt
import spatial_stats
import time
import tempfile
from linear_powerspect import get_powerspectra
import pickle as pkl
import sys

# This shows how to use the code to run a Zeldovich simulation with the power spectrum in the directory
# It takes a redshift as an input parameter, and an optional smoothing window parameter (default = 1 Mpc/h)
# It uses a boxsize of 1Gpc/h and 128^3 particles and redshift-space distortion parameter f=0.5
# computes the 1d and 2d correlation functions and plots both

def test(redshift, smw=1.0):
    pkinit = N.loadtxt('pk_indra7313.txt')
    boxsize = 700.0
    ngrid = 128

    dens0, x, y, z = Z.run(redshift, pkinit, 0.5, boxsize=boxsize,
                           ngrid=ngrid, exactpk=True, smw=smw, seed=N.random.randint(1000))
    dens = cic_dens_wrapper.get_dens(x, y, z, ngrid, boxsize)

    k, pk = spatial_stats.getPk(
        dens, nkbins=40, boxsize=boxsize, deconvolve_cic=True, exp_smooth=0.0)
    plt.figure(1)
    plt.loglog(k, pk, marker='*')
    # plt.xlim([0, 200])


def get_zel_powerspects(redshift, om_m, om_b, h, sigma_8, ns, n_runs, boxsize=350.0, ngrid=128, n_bins=40, smw=1.0, verbose=True, return_linear=False, save=None):
    pkinit = get_powerspectra(om_m, om_b, h, sigma_8, ns)

    ks, pks = [], []
    for i in xrange(n_runs):
        dens0, x, y, z = Z.run(redshift, pkinit, 0.5, boxsize=boxsize,
                               ngrid=ngrid, exactpk=False, smw=smw, seed=N.random.randint(1000000))
        dens = cic_dens_wrapper.get_dens(x, y, z, ngrid, boxsize)

        k, pk = spatial_stats.getPk(
            dens, nkbins=n_bins, boxsize=boxsize, deconvolve_cic=True, exp_smooth=0.0)
        ks.append(k)
        pks.append(pk)
        if verbose and i!=n_runs-1:
            print "Done {} out of {} runs".format(i, n_runs)
    print N.shape(pks)
    if save is not None:
        pkl.dump({'om_m':om_m, 'om_b':om_b, 'h':h, 'sigma_8':sigma_8, 'ns':ns, 'ks':N.array(ks).transpose(), 'pks':N.array(pks).transpose()}, open(save,'w'))
    if return_linear is False:
        return N.array(ks).transpose(), N.array(pks).transpose()
    else:
        return N.array(ks).transpose(), N.array(pks).transpose(), pkinit

def get_one_powerspectrum(om_m, om_b, h, sigma_8, ns, redshift, boxsize, ngrid, n_bins, smw=1.0):
    pkinit = get_powerspectra(om_m, om_b, h, sigma_8, ns)
    dens0, x, y, z = Z.run(redshift, pkinit, 0.5, boxsize=boxsize,
                           ngrid=ngrid, exactpk=False, smw=smw, seed=N.random.randint(1000000))
    dens = cic_dens_wrapper.get_dens(x, y, z, ngrid, boxsize)

    k, pk = spatial_stats.getPk(
        dens, nkbins=n_bins, boxsize=boxsize, deconvolve_cic=True, exp_smooth=0.0)

    return k, pk


def make_two_d_parameter_space(number_of_runs, filename, sig8_low=0.45, sig8_high=1.05, om_m_low=0.2, om_m_high=0.7, om_b=0.04, h=0.7, ns=1.0, bins=40, verbose=True):
    redshift=0.0
    boxsize=350.0
    ngrid=128
    smw=1.0
    input_dict = {}
    output_dict = {}
    rest_of_cosm = {'h':h, 'om_b':om_b, 'ns':ns}
    for i in xrange(number_of_runs):
        om_m_thisone= (om_m_high-om_m_low)*N.random.rand() + om_m_low
        sig_8_thisone= (sig8_high - sig8_low) * N.random.rand() + sig8_low
        input_dict[i] = [sig_8_thisone, om_m_thisone]
        output_dict[i] = get_one_powerspectrum(om_m_thisone, om_b, h, sig_8_thisone, ns, redshift, boxsize, ngrid, bins)
        if verbose and i!=number_of_runs-1:
            print "Done {} out of {} runs".format(i, number_of_runs)
    pkl.dump((input_dict, output_dict, rest_of_cosm), open(filename, 'w'))


if __name__ == '__main__':
    if sys.argv[1] == 'production':
        bins = int(sys.argv[2])
        make_two_d_parameter_space(10000, 'data/{}_bins_run.pkl'.format(bins), bins=bins)

    elif sys.argv[1] == 'nice_run':
        iterboi = int(sys.argv[2])
        n_bins = int(sys.argv[3])
        sig8_low=0.45
        sig8_high=1.05
        om_m_low=0.2
        om_m_high=0.7
        om_b=0.04
        h=0.7
        ns=1.0
        redshift=0.
        n_runs=10000
        om_m= (om_m_high-om_m_low)*N.random.rand() + om_m_low
        sigma_8 = (sig8_high - sig8_low) * N.random.rand() + sig8_low
        k2, pk2, pkinit2 = get_zel_powerspects(redshift, om_m, om_b, h, sigma_8, ns, n_runs, return_linear=True, save='data/full_run_number={}_bins={}.pkl'.format(iterboi, n_bins), n_bins=n_bins)
    elif sys.argv[1] == 'test':
        om_m = 0.25
        om_b = 0.06
        h = 0.6
        sigma_8 = 0.8
        ns =0.9
        redshift=0.0
        boxsize=350.0
        ngrid = 128
        k, pk = get_one_powerspectrum(om_m, om_b, h, sigma_8, ns, redshift, boxsize, ngrid, 10, smw=1.0)
        plt.loglog(k, pk)
        k, pk = get_one_powerspectrum(om_m, om_b, h, sigma_8, ns, redshift, boxsize, ngrid, 100, smw=1.0)
        plt.loglog(k, pk)
        plt.show()
    # make_two_d_parameter_space(10, 'data/first_test.pkl')

    # k2, pk2, pkinit2 = get_zel_powerspects(0.0, 0.25, 0.04, 0.7, 0.9, 1.0, 10000, return_linear=True, save='data/high_sigma_8_large.pkl')

    # k, pk, pkinit = get_zel_powerspects(0.0, 0.25, 0.04, 0.7, 0.6, 1.0, 10000, return_linear=True, save='data/low_sigma_8_large.pkl')
    # make_two_d_parameter_space(10000, 'data/many_runs.pkl')

    # print k.shape,pk.shape
    # plt.loglog(k, pk)
    # plt.loglog(k2, pk2)
    # # plt.loglog(pkinit[:,0], pkinit[:,1])
    # # plt.loglog(pkinit2[:,0], pkinit2[:,1])
    # # plt.figure()
    # # plt.imshow(N.corrcoef(pk),vmax=1.0, vmin=-1.0)
    # # plt.colorbar()
    # print pkl.load(open('data/test.pkl','r'))
    # plt.show()