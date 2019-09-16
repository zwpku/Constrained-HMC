#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

def U(theta):
  x0 = c * math.cos(theta)
  return 0.0
#  return 2.0 * x0**2

dim = 3
output_every_k = 1
R = 1.0
r = 0.5
lc = ['k', 'r', 'k', 'c', 'm', 'y', 'y']
ls = ['--', '--', '--', '--', '--', '--', '--']
mks= ['o', 's', 'h', 'x', 'v', 'D', '^']
label_names = ['PR', 'PR-far', 'PR50', 'Newton', 'Hom', 'Hom-far', 'Hom50']
job_idx_vec = [2, 6, 4, 1, 3, 7, 5]
mk_every_vec = [10, 12]
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

# compute the reference density 
num_x = 500
dtheta = 2 * math.pi / num_x 
true_density = [0] * num_x
theta_vec = [(i + 0.5) * dtheta for i in range(num_x)]
norm_z = 0
for i in range(num_x):
  true_density[i] = 1.0 / (2 * math.pi) * (1.0 + r * 1.0 / R * math.cos(theta_vec[i]))
  norm_z += true_density[i] * dtheta

#normalize the density
true_density = [true_density[i] / norm_z for i in range(num_x)]

# plot reference density
plt.plot(theta_vec, true_density, linestyle='-', color='k', label='true')

for idx in range(7):
    job_id = job_idx_vec[idx]
    working_dir_name = '../working_dir_task%d/' % (job_id)

    data_file_name = '%s/qoi_counter_%d.txt' % (working_dir_name, job_id)
    infile = open(data_file_name, 'r')

    N, num_qoi = [int(x) for x in infile.readline().split()]
    qoi_hist_info = [[float(x) for x in infile.readline().split()] for idx in range(num_qoi)]
    qoi_counter = [[int(float(x)) for x in infile.readline().split()] for idx in range(num_qoi)]

    print ("%d QoI are loaded from: %s" % (num_qoi, data_file_name))
    xx = np.linspace( qoi_hist_info[0][1], qoi_hist_info[0][2], int(qoi_hist_info[0][0]) )
    density = [x  / (N * qoi_hist_info[0][3]) for x in qoi_counter[0]]

    plt.plot(xx, density, color=lc[0], linestyle=ls[idx], marker=mks[idx],
                markevery= 10, fillstyle='none', label= label_names[idx] )


plt.xticks([0, math.pi/2, math.pi, math.pi * 3.0 /2, 2.0 * math.pi], ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'], fontsize=22)
plt.yticks(fontsize=22)
plt.gca().set_ylim(0.03, 0.34)
plt.legend(bbox_to_anchor=(0.17, 0.5, 0.9, 0.52), fontsize=20, ncol=2)
fig.tight_layout()
out_fig_name = './hist_angle_3dtorus.eps' 
fig.savefig(out_fig_name)

