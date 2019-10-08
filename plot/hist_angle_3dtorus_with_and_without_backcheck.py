#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

def V(phi, theta):
  return 0.0

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

job_id_vec = [1, 0, 0, 0]
working_dir_vec = ['../working_dir_task1/', '../working_dir_task0-2/', '../working_dir_task0/', '../working_dir_task0-1/']

dim = 3
output_every_k = 1
R = 1.0
r = 0.5
beta = 20.0

lc = ['k', 'r', 'k', 'c', 'm', 'y', 'y']
ls = ['--', '-.', ':', '--']
mks= ['+', 's', 'D', 'o', 'x', 'v', '^']
label_names = [r'backward, $\tau=0.5$', r'backward, $\tau=1.0$', r'no backward, $\tau=0.5$', r'no backward, $\tau=1.0$']
lab_vec = [r'$\phi$', r'$\theta$']

for angle_idx in range(2):
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ref_density = reference_density(angle_idx)
    xx = np.linspace(0, 2.0 * math.pi, len(ref_density) )
    plt.plot(xx, ref_density, linestyle='-', color='k', label='true')

#    for j in range(2, len(working_dir_vec)):
    for j in range(0, 2):
        data_file_name = '%s/qoi_counter_%d.txt' % (working_dir_vec[j], job_id_vec[j])
        infile = open(data_file_name, 'r')

        N, num_qoi = [int(x) for x in infile.readline().split()]
        qoi_hist_info = [[float(x) for x in infile.readline().split()] for idx in range(num_qoi)]
        qoi_counter = [[int(float(x)) for x in infile.readline().split()] for idx in range(num_qoi)]

        print ("%d QoI are loaded from: %s" % (num_qoi, data_file_name))
        xx = np.linspace( qoi_hist_info[0][1], qoi_hist_info[0][2], int(qoi_hist_info[0][0]) )
        density = [x  / (N * qoi_hist_info[0][3]) for x in qoi_counter[angle_idx]]

        plt.plot( xx, density, color=lc[0], linestyle=ls[j], marker=mks[j], markevery= 20, fillstyle='none', label= label_names[j] )

    plt.xticks([0, math.pi/2, math.pi, math.pi * 3.0 /2, 2.0 * math.pi], ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'], fontsize=22)
    plt.yticks(fontsize=22)
    plt.title('Probability density of '+lab_vec[angle_idx], fontsize=22)

    plt.gca().set_ylim(0.05, 0.34)
    plt.legend(bbox_to_anchor=(0.30, 0.75, 0.4, 0.22), fontsize=19, ncol=1, frameon=False)
    fig.tight_layout()
    out_fig_name = './cmp_back_check_%dth.eps' % angle_idx
    fig.savefig(out_fig_name)

