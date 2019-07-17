#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

def U(theta):
  x0 = c * math.cos(theta)
  return 2.0 * x0**2

working_dir_name = '../'
job_id = 1
dim = 3
output_every_k = 1
c = 3.0
c2 = c**2
lc = ['b', 'r', 'k', 'c', 'm', 'y']

data_file_name = '%s/data/data_%d.txt' % (working_dir_name, job_id)

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
xv_data = np.loadtxt(data_file_name)
N = len(xv_data[:,0])
print ("%d samples are loaded from: %s" % (N, data_file_name))

angle_vec_from_data = [math.atan2(xv_data[i][1], xv_data[i][0]/c) for i in range(N)] 
for i in range(N):
  if angle_vec_from_data[i] < 0:
      angle_vec_from_data[i] += 2 * math.pi

plt.plot(angle_vec_from_data[::output_every_k], color=lc[0], linestyle='-', label=r'$\theta$')
fig.tight_layout()
ax.legend(bbox_to_anchor=(0.5, 0, 0.5, 0.5))
out_fig_name = '%s/fig/traj_ellipse_angle_%d.eps' % (working_dir_name, job_id)
fig.savefig(out_fig_name)

plt.clf()

num_x = 500
dtheta = 2 * math.pi / num_x 
true_density = [0] * num_x
theta_vec = [(i + 0.5) * dtheta for i in range(num_x)]
norm_z = 0
for i in range(num_x):
  true_density[i] = math.exp(-U(theta_vec[i])) * math.sqrt(c2 * math.sin(theta_vec[i])**2 + math.cos(theta_vec[i])**2)
  norm_z += true_density[i] * dtheta
true_density = [true_density[i] / norm_z for i in range(num_x)]
  
plt.hist(angle_vec_from_data, bins=500, density=True, histtype='step',color=lc[0], label='hist')
plt.plot(theta_vec, true_density,color=lc[1], linestyle='-', label='true density')
fig.tight_layout()
ax.legend(bbox_to_anchor=(0.5, 0, 0.5, 0.5))
out_fig_name = '%s/fig/hist_ellipse_angle_%d.eps' % (working_dir_name, job_id)
fig.savefig(out_fig_name)
