#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

job_id = 1
working_dir_vec = ['../working_dir_ex2_newton_test/', '../working_dir_ex2_homotopy_potential_test/', '../working_dir_ex2_homotopy10_potential_test/']
dim = 10
lc = ['k', 'r', 'k', 'c', 'm', 'y']
mks= [None, '^', 'x', 'o', 'x', 'v', '^']
ls = ['--', '-', '-', '--', '--', '--', '--']
num_of_qoi = 3
label_name = [r'Newton', r'Hom', r'Hom10']
qoi_index_vec = [1, 2, 4]
oft = [10, 10, 10] 

fig, ax = plt.subplots(1, 3, sharey=True, gridspec_kw={'hspace': 0}, figsize=(14, 4))

# loop for each scheme
for idx in range(3):
    working_dir_name = working_dir_vec[idx]
    data_file_name = '%s/qoi_counter_%d.txt' % (working_dir_name, job_id)
    infile = open(data_file_name, 'r')
    N, num_qoi = [int(x) for x in infile.readline().split()]
    qoi_hist_info = [[float(x) for x in infile.readline().split()] for idx in range(num_qoi)]
    qoi_counter = [[int(float(x)) for x in infile.readline().split()] for idx in range(num_qoi)]
    print ("%d QoI are loaded from: %s" % (num_qoi, data_file_name))

    for j in range(3):
        q_idx = qoi_index_vec[j]
        xx = np.linspace( qoi_hist_info[q_idx][1], qoi_hist_info[q_idx][2], int(qoi_hist_info[q_idx][0]) )
        density = [x  / (N * qoi_hist_info[q_idx][3]) for x in qoi_counter[q_idx]]
        ax[j].plot(xx, density, color=lc[0], linestyle=ls[idx], marker=mks[idx], markersize=8, markevery= oft[idx], fillstyle='none', label=label_name[idx])
        ax[j].set_xticks([-3, -2, -1, 0, 1,2,3])
        ax[j].set_xticklabels([-3, -2, -1, 0, 1,2,3], fontsize=18)
        ax[j].set_title(r'$x_%d$' % q_idx, fontsize=18)
#        ax[j].legend(bbox_to_anchor=(0.2, 0.6, 0.4, 0.2), fontsize=12, frameon=False)
        ax[j].legend(loc="upper left", fontsize=12, frameon=False)

ax[0].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
ax[0].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8], fontsize=22)
plt.gca().set_ylim(0.0, 0.9)
#fig.tight_layout()
out_fig_name = './ex2_hist_qoi.eps' 
fig.savefig(out_fig_name)

print ("Figure saved to : %s" % out_fig_name)

