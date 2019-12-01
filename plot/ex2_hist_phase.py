#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

job_id = 1
working_dir_name = '../working_dir_ex2_homotopy50_potential_test/' 
dim = 10
beta = 1.0
lc = ['b', 'r', 'k', 'c', 'm', 'y']
label_name = ['phase idx']

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

#    plt.xticks([0, math.pi/2, math.pi, math.pi * 3.0 /2, 2.0 * math.pi], ['0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'], fontsize=22)
    plt.yticks(fontsize=18)
    plt.gca().set_ylim(0.0, 0.4)
    plt.legend(bbox_to_anchor=(0.3, 0.5, 0.5, 0.5), fontsize=18)
    fig.tight_layout()
    out_fig_name = '%s/hist_qoi_%dth.eps' % (working_dir_name, idx)
    fig.savefig(out_fig_name)

