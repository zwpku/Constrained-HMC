#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

job_id = 1
N = int(1e7)
working_dir_vec = ['../working_dir_ex2_newton_test-new/', '../working_dir_ex2_homotopy_potential_test-new/', '../working_dir_ex2_homotopy10_potential_test-new/']
dim = 10
lc = ['k', 'r', 'k', 'c', 'm', 'y']
mks= [None, '^', 'x', 'o', 'x', 'v', '^']
ls = ['--', '-', '-', '--', '--', '--', '--']
scheme_name = [r'Newton', r'Hom', r'Hom10']
oft = [1, 5000, 5000]

fig, ax = plt.subplots(1, 3, sharey=True, gridspec_kw={'hspace': 0}, figsize=(14, 4))

for idx in range(3):
    working_dir_name = working_dir_vec[idx]
    qoi_data_file_name = '%s/qoi_data_%d.txt' % (working_dir_name, job_id)
    print("loading phase data from : ", qoi_data_file_name)
    xv_data = np.loadtxt(qoi_data_file_name, usecols=0)
    len_data = len(xv_data)
    print ("%d samples are loaded from: %s" % (len_data, qoi_data_file_name))
    output_every_k = oft[idx]
    num_data_to_plot = int(len_data / output_every_k)
    xx = np.linspace(0, N, num_data_to_plot)
    ax[idx].plot(xx, xv_data[0::output_every_k], '-', color='b', linewidth=1)
    ax[idx].set_ylim([-0.5, 3.5])
    ax[idx].set_yticks([0, 1, 2, 3])
    ax[idx].set_yticklabels([r'$\mathcal{C}_0$',r'$\mathcal{C}_1$',r'$\mathcal{C}_2$',r'$\mathcal{C}_3$'], fontsize=22)
    ax[idx].set_xticks([0, N/2, N])
    ax[idx].set_xticklabels([r'$0$', r'$5\times 10^6$', r'$10^7$'], fontsize=18)
    ax[idx].set_title(scheme_name[idx], fontsize=18)

out_fig_name = './ex2_traj_phase.eps' 
fig.savefig(out_fig_name)

