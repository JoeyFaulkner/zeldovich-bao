import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np

a = pkl.load(open('data/full_run_number=1_bins=40.pkl', 'r'))
print a['pks']

# plt.loglog(a['ks'][:,0], a['pks'][:, 0])
# print a['pks'][:, :1000]/
true = np.corrcoef(a['pks'][:, :])
for bins in [500, 1000, 2000, 3000, 5000, 10000]:
    plt.imshow(true - np.corrcoef(a['pks'][:, :bins]), vmax=0.1, vmin=-0.1)
    print np.mean(np.abs(true - np.corrcoef(a['pks'][:, :bins])))
    plt.colorbar()
    plt.show()