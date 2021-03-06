#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

#working_dir_name = '../working_dir_ex2_homotopy_test'
working_dir_name = '../working_dir_ex2_homotopy_potential_test'
#working_dir_name = '../working_dir_ex2_homotopy50_potential_test'
job_id = 1
dim = 10
lc = ['b', 'r', 'k', 'c', 'm', 'y']

data_file_name = '%s/traj_data_%d.txt' % (working_dir_name, job_id)

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
xv_data = np.loadtxt(data_file_name)
print ("%d samples are loaded from: %s" % (len(xv_data[:,0]), data_file_name))

# plot the position
for i in range(dim):
    plt.clf()
    plt.hist(xv_data[:, i], bins=100, density=True)
    fig.tight_layout()
    out_fig_name = '%s/hist_plot_x%d_%d.eps' % (working_dir_name, i, job_id)
    fig.savefig(out_fig_name)

# plot the velocity
for i in range(dim):
    plt.clf()
    plt.hist(xv_data[:,dim + i], bins=100, density=True)
#    ax.legend(bbox_to_anchor=(0.5, 0, 0.5, 0.5))
    fig.tight_layout()
    out_fig_name = '%s/hist_plot_v%d_%d.eps' % (working_dir_name, i, job_id)
    fig.savefig(out_fig_name)

