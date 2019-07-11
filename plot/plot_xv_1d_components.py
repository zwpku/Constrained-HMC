#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

working_dir_name = '../'
dim = 2
output_every_k = 3
lc = ['b', 'r', 'k', 'c', 'm', 'y']

data_file_name = '%s/data/data.txt' % (working_dir_name)

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
xv_data = np.loadtxt(data_file_name)
print ("%d samples are loaded from: %s" % (len(xv_data[:,0]), data_file_name))

# plot the position
for i in range(dim):
    plt.clf()
    plt.plot(xv_data[::output_every_k,i], color=lc[i], linestyle='-', label='x_%d' % i)
#    ax.legend(bbox_to_anchor=(0.5, 0, 0.5, 0.5))
#    ax.set_ylim(-3.0, 3.00)
#    ax.set_yticks(np.arange(-3.0, 3.0, step=1.0))
    fig.tight_layout()
    out_fig_name = '%s/fig/traj_plot_x%d.eps' % (working_dir_name, i)
    fig.savefig(out_fig_name)

# plot the velocity
for i in range(dim):
    plt.clf()
    plt.plot(xv_data[::output_every_k,dim + i], color=lc[i], linestyle='-', label='v_%d' % i)
    fig.tight_layout()
    out_fig_name = '%s/fig/traj_plot_v%d.eps' % (working_dir_name, i)
    fig.savefig(out_fig_name)
