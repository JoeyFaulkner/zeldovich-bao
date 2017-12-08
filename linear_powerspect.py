import camb
import numpy as np


# def write_powerspectrum(om_m, om_b, h, sigma_8, ns, filename):
#     k, pk = get_powerspectra(om_m, om_b, h, sigma_8, ns)
#     out_bit = np.vstack((k,pk)).transpose()
#     np.savetxt(filename, out_bit)


def get_powerspectra(om_m, om_b, h, sigma_8, ns):
    om_c = om_m - om_b
    initial = camb.initialpower.InitialPowerParams()
    initial.set_params(ns=ns)
    model = camb.model.CAMBparams()
    model.set_initial_power(initial)
    model.set_cosmology(H0=100 * h, ombh2=om_b * h**2, omch2=om_c * h**2)
    model.set_matter_power(k_per_logint=0, redshifts=[0])
    model.validate()

    data = camb.CAMBdata()
    data.set_params(model)

    q = data.get_matter_power_spectrum(
        params=model, minkh=0.0001, maxkh=32.126, npoints=3000)
    sig_8out = data.get_matter_transfer_data().sigma_8[0][0]

    k = q[0]
    pk = q[2][0]



    pk_out = pk

    return np.vstack((k, pk_out * (sigma_8 ** 2 / sig_8out ** 2))).transpose()


if __name__ == '__main__':
    import tempfile
    f=tempfile.NamedTemporaryFile()
    import matplotlib.pyplot as plt
    write_powerspectrum(0.272, 0.0455, 0.702, 0.807, 0.961, f.name)

    ks, pk =np.loadtxt(f.name).transpose()

    plt.loglog(ks, pk )
    plt.show()