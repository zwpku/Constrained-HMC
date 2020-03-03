#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

def V(phi, theta):
  if pot_id == 0 :
      return 0.0
  else :
      x = [(R+r*math.cos(phi))*math.cos(theta), (R+r*math.cos(phi))*math.sin(theta), r * math.sin(phi)]
      return (x[0] - x[1])**2 + 5.0 * ((x[0]**2 + x[1]**2)/(R+r)**2 - 1)**2

def reference_density(idx) :
    # compute the reference density 
    num_x = 500
    dx = 2 * math.pi / num_x 
    true_density = [0] * num_x
    angle_vec = [(i + 0.5) * dx for i in range(num_x)]
    norm_z = 0
    for i in range(num_x):
      sum_tmp = 0.0
      # compute density for phi
      if idx == 0 :
          for j in range(num_x):
              arg = (j + 0.5) * dx
              sum_tmp += math.exp(-beta * V(angle_vec[i], arg)) * dx
          true_density[i] = (1.0 + r * 1.0 / R * math.cos(angle_vec[i])) * sum_tmp
      else :
      # compute density for theta
          for j in range(num_x):
              arg = (j + 0.5) * dx
              sum_tmp += (1.0 + r * 1.0 / R * math.cos(arg)) * math.exp(-beta * V(arg, angle_vec[i])) * dx
          true_density[i] = sum_tmp

      norm_z += true_density[i] * dx

    #normalize the density
    true_density = [true_density[i] / norm_z for i in range(num_x)]
    return true_density

job_id = 10
#working_dir_name = '../working_dir_task%d-dt1.0/' % job_id
working_dir_name = '../working_dir_task%d/' % job_id
#working_dir_name = '../working_dir_ex1-alpha0.5-task%d/' % job_id
dim = 3
pot_id = 1
output_every_k = 1
R = 1.0
r = 0.5
beta = 20.0
lc = ['b', 'r', 'k', 'c', 'm', 'y']
label_name = [r'$\phi$', r'$\theta$']

data_file_name = '%s/qoi_counter_%d.txt' % (working_dir_name, job_id)
infile = open(data_file_name, 'r')

N, num_qoi = [int(x) for x in infile.readline().split()]

qoi_hist_info = [[float(x) for x in infile.readline().split()] for idx in range(num_qoi)]
qoi_counter = [[int(float(x)) for x in infile.readline().split()] for idx in range(num_qoi)]
print ("%d QoI are loaded from: %s" % (num_qoi, data_file_name))

for idx in range(num_qoi):
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    xx = np.linspace( qoi_hist_info[0][1], qoi_hist_info[0][2], int(qoi_hist_info[0][0]) )
    density = [x  / (N * qoi_hist_info[0][3]) for x in qoi_counter[idx]]
    plt.plot(xx, density, color=lc[0], linestyle='-', label=label_name[idx])

    ref_density = reference_density(idx)
    xx = np.linspace( qoi_hist_info[0][1], qoi_hist_info[0][2], len(ref_density) )
    plt.plot(xx, ref_density, color=lc[0], linestyle='--', label='true')

    plt.xticks([0, math.pi/2, math.pi, math.pi * 3.0 /2, 2.0 * math.pi], ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'], fontsize=22)
    plt.yticks(fontsize=18)
#    plt.gca().set_ylim(0.0, 4.0)
    plt.legend(bbox_to_anchor=(0.3, 0.5, 0.5, 0.5), fontsize=18)
    fig.tight_layout()
    out_fig_name = '%s/hist_3dtorus_job%d_%dth.eps' % (working_dir_name, job_id, idx)
    fig.savefig(out_fig_name)

