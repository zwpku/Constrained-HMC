#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

job_id = 1
N = int(1e7)
working_dir_vec = ['../working_dir_ex2_newton_test-new/', '../working_dir_ex2_homotopy_potential_test-new/', '../working_dir_ex2_homotopy10_potential_test-new/']
dim = 10
lc = ['k', 'r', 'k', 'c', 'm', 'y']
mks= [None, '^', 'x', 'o', 'x', 'v', '^']
ls = ['--', '-', '-', '--', '--', '--', '--']
scheme_name = [r'Newton', r'Hom', r'Hom10']
oft = [1, 5000, 5000]

fig, ax = plt.subplots(1, 3, sharey=True, gridspec_kw={'hspace': 0}, figsize=(14, 4))

for idx in range(3):
    working_dir_name = working_dir_vec[idx]
    qoi_data_file_name = '%s/qoi_data_%d.txt' % (working_dir_name, job_id)
    print("loading phase data from : ", qoi_data_file_name)
    xv_data = np.loadtxt(qoi_data_file_name, usecols=0)
    len_data = len(xv_data)
    print ("%d samples are loaded from: %s" % (len_data, qoi_data_file_name))
    ind_vec = np.where(np.fabs(xv_data[:-1] - xv_data[1:])>0.5)[0] 
    num_pc = len(ind_vec)
    print("num of phase changes: %d" % num_pc)

    counter = np.zeros((4,4))
    for i in range(num_pc):
        ii = int(xv_data[ind_vec[i]])
        jj = int(xv_data[1 + ind_vec[i]])
        counter[ii][jj] += 1

    for i in range(4):
        for j in range(4):
            print("%.2e  " % (counter[i][j] * 1.0 / N), end="")
        print(end="\n")


