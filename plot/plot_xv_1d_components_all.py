#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

dim = 3
output_every_k = [1, 1800, 1800]
lc = ['k', 'r', 'k', 'c', 'm', 'y']
ls = ['--', '--', '--', '--', '--', '--', '--']
mks= ['+', 's', 'o', 'x', 'v', 'D', 'h']
lab_vec = [r'$\phi$', r'$\theta$']

title_names = ['Newton', 'PR', 'PR50-far']
job_idx_vec = [8, 9, 11]
angle_idx = 1
N = int(1e7)
num_jobs = len(job_idx_vec)

fig, ax = plt.subplots(1, num_jobs, sharey=True, gridspec_kw={'hspace': 0}, figsize=(8, 3))

for idx in range(num_jobs):
    job_id = job_idx_vec[idx]
    working_dir_name = '../working_dir_task%d/' % job_id
    qoi_data_file_name = '%s/qoi_data_%d.txt' % (working_dir_name, job_id)
    xv_data = np.loadtxt(qoi_data_file_name)
    len_data = len(xv_data[:,0])
    print ("%d samples are loaded from: %s" % (len_data, qoi_data_file_name))
    num_data_to_plot = int(len_data / output_every_k[idx])
    xx = np.linspace(0, N, num_data_to_plot)
    ax[idx].plot(xx, xv_data[0::output_every_k[idx], angle_idx][0:num_data_to_plot], color='b', linestyle='-', linewidth=1)
    ax[idx].set_ylim([0, 3.0 / 2 * math.pi])
    ax[idx].set_xticks([0, N/2, N])
    ax[idx].set_xticklabels([r'$0$', r'$10^6$', r'$10^7$'], fontsize=18)
    ax[idx].set_title(title_names[idx], fontsize=18)

#    plt.plot(xv_data[1::,i], color=lc[i], linestyle='-', label='x_%d' % i)
#    ax.legend(bbox_to_anchor=(0.5, 0, 0.5, 0.5))
#    ax.set_ylim(-3.0, 3.00)
#    ax.set_yticks(np.arange(-3.0, 3.0, step=1.0))

ax[0].set_yticks([math.pi/4, math.pi * 3.0 / 4, math.pi * 5.0 /4])
ax[0].set_yticklabels([r'$\frac{\pi}{4}$', r'$\frac{3\pi}{4}$', r'$\frac{5\pi}{4}$'], fontsize=22)
ax[0].set_ylabel(r'$\theta$', fontsize=18)

#fig.tight_layout()
out_fig_name = './traj_plot_angle_%dth.eps' % (angle_idx)
fig.savefig(out_fig_name)

