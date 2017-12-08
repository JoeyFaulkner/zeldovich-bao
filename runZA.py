import zeldovich as Z
import cic_dens_wrapper
import numpy as N
import matplotlib.pyplot as plt
import spatial_stats
import time
import tempfile
from linear_powerspect import get_powerspectra
import pickle as pkl

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
                               ngrid=ngrid, exactpk=True, smw=smw, seed=N.random.randint(1000000))
        dens = cic_dens_wrapper.get_dens(x, y, z, ngrid, boxsize)

        k, pk = spatial_stats.getPk(
            dens, nkbins=n_bins, boxsize=boxsize, deconvolve_cic=True, exp_smooth=0.0)
        ks.append(k)
        pks.append(pk)
        if verbose and i!=n_runs-1:
            print "Done {} out of {} runs".format(i, n_runs)

    if save is not False:
        pkl.dump((N.array(ks).transpose(), N.array(pks).transpose()), open(save,'w'))
    if return_linear is False:
        return N.array(ks).transpose(), N.array(pks).transpose()
    else:
        return N.array(ks).transpose(), N.array(pks).transpose(), pkinit
if __name__ == '__main__':
    k, pk, pkinit = get_zel_powerspects(0.0, 0.30, 0.04, 0.7, 0.9, 1.0, 1, return_linear=True)
    k2, pk2, pkinit2 = get_zel_powerspects(0.0, 0.25, 0.04, 0.7, 0.9, 1.0, 1, return_linear=True)

    print k.shape,pk.shape
    plt.loglog(k, pk)
    plt.loglog(k2, pk2)
    # plt.loglog(pkinit[:,0], pkinit[:,1])
    # plt.loglog(pkinit2[:,0], pkinit2[:,1])
    # plt.figure()
    # plt.imshow(N.corrcoef(pk),vmax=1.0, vmin=-1.0)
    # plt.colorbar()
    plt.show()