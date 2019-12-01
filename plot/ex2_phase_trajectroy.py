#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

dim = 10
output_every_k = 1
lc = ['k', 'r', 'k', 'c', 'm', 'y']
ls = ['--', '--', '--', '--', '--', '--', '--']
mks= ['+', 's', 'o', 'x', 'v', 'D', 'h']

N = int(1e5)

fig, ax = plt.subplots(1, 1, sharey=True, gridspec_kw={'hspace': 0}, figsize=(10, 3))

job_id = 1
#working_dir_name = '../working_dir_ex2_test/' 
working_dir_name = '../working_dir_ex2_homotopy_potential_test/' 
qoi_data_file_name = '%s/qoi_data_%d.txt' % (working_dir_name, job_id)
xv_data = np.loadtxt(qoi_data_file_name)
len_data = len(xv_data)
print ("%d samples are loaded from: %s" % (len_data, qoi_data_file_name))
num_data_to_plot = int(len_data / output_every_k)
print (num_data_to_plot)
xx = np.linspace(0, N, num_data_to_plot)
ax.plot(xx, xv_data[0::output_every_k][0:num_data_to_plot], color='b', linestyle='-', linewidth=1)
ax.set_ylim([-1, 4.0])

out_fig_name = '%s/traj_phase.eps' % working_dir_name
fig.savefig(out_fig_name)

