#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

lc = ['b', 'r', 'k', 'c', 'm', 'y']

fig, ax = plt.subplots(1, 1, figsize=(8, 5))

t = np.array(np.linspace(0, math.pi * 2, 100))
plt.plot(np.sqrt(4 + np.sin(t)), np.sqrt(4 + np.cos(t)), color=lc[0], linestyle='-')
plt.plot(-np.sqrt(4 + np.sin(t)), np.sqrt(4 + np.cos(t)), color=lc[0], linestyle='-')
plt.plot(np.sqrt(4 + np.sin(t)), -np.sqrt(4 + np.cos(t)), color=lc[0], linestyle='-')
plt.plot(-np.sqrt(4 + np.sin(t)), -np.sqrt(4 + np.cos(t)), color=lc[0], linestyle='-')
#    ax.legend(bbox_to_anchor=(0.5, 0, 0.5, 0.5))
ax.set_xlim(-3.0, 3.00)
ax.set_ylim(-3.0, 3.00)
#    ax.set_yticks(np.arange(-3.0, 3.0, step=1.0))
fig.tight_layout()
out_fig_name = './test.eps' 
fig.savefig(out_fig_name)

