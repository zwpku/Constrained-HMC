#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

working_dir_name = '../'
job_id = 3
dim = 2
output_every_k = 1
lc = ['b', 'r', 'k', 'c', 'm', 'y']

data_file_name = '%s/data/data_%d.txt' % (working_dir_name, job_id)

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
xv_data = np.loadtxt(data_file_name)
N = len(xv_data[:,0])
print ("%d samples are loaded from: %s" % (N, data_file_name))

distance_vec = [np.linalg.norm(xv_data[i,0:dim] - xv_data[i-1,0:dim]) for i in range(1, N)] 

xx = range(1, 5001, output_every_k)
plt.plot(xx, distance_vec[0:5000][::output_every_k], color=lc[0], linestyle='-', label=r'distance')

fig.tight_layout()
out_fig_name = '%s/fig/traj_jump_distance_%d.eps' % (working_dir_name, job_id)
fig.savefig(out_fig_name)

plt.clf()

histdata, histedge = np.histogram(distance_vec, bins=100, density=True)
histcenter = (histedge[:-1] + histedge[1:]) * 0.5

plt.plot(histcenter, histdata, color=lc[0], linestyle='-', label='distance')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(bbox_to_anchor=(0.3, 0.5, 0.5, 0.5), fontsize=18)
fig.tight_layout()
out_fig_name = '%s/fig/hist_jump_distance_%d.eps' % (working_dir_name, job_id)
fig.savefig(out_fig_name)
