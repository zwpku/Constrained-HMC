#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

def U(theta):
  x0 = c * math.cos(theta)
  return 0.0
#  return 2.0 * x0**2

working_dir_name = '../working_dir_task4/'
job_id = 4
dim = 3
output_every_k = 1
R = 1.0
r = 0.5
lc = ['b', 'r', 'k', 'c', 'm', 'y']

data_file_name = '%s/qoi_counter_%d.txt' % (working_dir_name, job_id)
infile = open(data_file_name, 'r')

fig, ax = plt.subplots(1, 1, figsize=(8, 5))

N, num_qoi = [int(x) for x in infile.readline().split()]

qoi_hist_info = [[float(x) for x in infile.readline().split()] for idx in range(num_qoi)]
qoi_counter = [[int(float(x)) for x in infile.readline().split()] for idx in range(num_qoi)]

print ("%d QoI are loaded from: %s" % (num_qoi, data_file_name))

xx = np.linspace( qoi_hist_info[0][1], qoi_hist_info[0][2], int(qoi_hist_info[0][0]) )
density = [x  / (N * qoi_hist_info[0][3]) for x in qoi_counter[0]]

plt.plot(xx, density, color=lc[0], linestyle='-', label=r'$\theta$')

num_x = 500
dtheta = 2 * math.pi / num_x 
true_density = [0] * num_x
theta_vec = [(i + 0.5) * dtheta for i in range(num_x)]
norm_z = 0
for i in range(num_x):
  true_density[i] = 1.0 / (2 * math.pi) * (1.0 + r * 1.0 / R * math.cos(theta_vec[i]))
  norm_z += true_density[i] * dtheta
true_density = [true_density[i] / norm_z for i in range(num_x)]

plt.plot(theta_vec, true_density, color=lc[1], linestyle='-', label='true density')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.gca().set_ylim(0.0, 0.3)
plt.legend(bbox_to_anchor=(0.3, 0.5, 0.5, 0.5), fontsize=18)
fig.tight_layout()
out_fig_name = '%s/hist_angle_3dtorus_%d.eps' % (working_dir_name, job_id)
fig.savefig(out_fig_name)
