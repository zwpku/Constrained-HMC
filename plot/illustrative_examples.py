#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import interpolate
from scipy.interpolate import spline

lc = ['k', 'r', 'k', 'c', 'm', 'y']
fs = 18

fig, ax = plt.subplots(1, 3, figsize=(12, 4))

t = np.linspace(0, np.pi*2, 100)
ax[0].plot(np.cos(t), np.sin(t), color=lc[0], linestyle='-')

ax[0].axhline(y=1,xmin=0.02,xmax=0.98, color=lc[0],linestyle='--', linewidth=1)
ax[0].axvline(x=-1.4,ymin=0.15,ymax=0.90, color=lc[0],linestyle='--',linewidth=1)
ax[0].axvline(x=-1,ymin=0.15,ymax=0.90, color=lc[0],linestyle='--',linewidth=1)
ax[0].axvline(x=-0.6,ymin=0.15,ymax=0.90, color=lc[0],linestyle='--',linewidth=1)
ax[0].axvline(x=1,ymin=0.15,ymax=0.90, color=lc[0],linestyle='--',linewidth=1)

ax[0].plot(-1.4,1, 'k.')
ax[0].plot(-1,1, 'k.')
ax[0].plot(-0.6,1, 'k.')
ax[0].plot(1,1, 'k.')
ax[0].plot(0,1, 'k.')

ax[0].plot(1,0, 'ko')
ax[0].plot(-1,0, 'ko')
ax[0].plot(-0.6, np.sqrt(1-0.6**2), 'ko')
ax[0].plot(-0.6, -np.sqrt(1-0.6**2), 'ko')

ax[0].arrow(0, 1, 0, 0.3, facecolor='black', width=0.01,head_width=0.05, head_length=0.1)
ax[0].text(0.1, 1.4, r'$\nabla\xi(x)$', fontsize=fs)
ax[0].annotate(r'$T_x\Sigma$', xy=(1.4, 1.1), xytext=(1.5, 1.6),
        arrowprops=dict(facecolor='black', shrink=0.05,width=0.6,headwidth=8.0,
            headlength=6), fontsize=fs)
ax[0].text(0.83,-0.0, r'$y$', fontsize=fs)
ax[0].text(1.25, -0.3, r'$\nabla\xi(y)$', fontsize=fs)
ax[0].arrow(1, 0, 0.3, 0.0, facecolor='black', width=0.01,head_width=0.05, head_length=0.1)
ax[0].text(0.00, 0.8, r'$x$', fontsize=fs)
ax[0].text(-1.35, 1.1, r'$v_1$', fontsize=fs)
ax[0].text(-0.95, 1.1, r'$v_2$', fontsize=fs)
ax[0].text(-0.55, 1.1, r'$v_3$', fontsize=fs)
ax[0].text(0.75, 1.1, r'$v_4$', fontsize=fs)
#ax[0].spines['left'].set_position('center')
ax[0].spines['left'].set_color('none')
ax[0].spines['bottom'].set_color('none')
# Eliminate upper and right axes
ax[0].spines['right'].set_color('none')
ax[0].spines['top'].set_color('none')
# Show ticks in the left and lower axes only
ax[0].xaxis.set_ticks_position('bottom')
ax[0].yaxis.set_ticks_position('left')
#plt.plot(histcenter, histdata, color=lc[0], linestyle='-', label='distance')
plt.xticks(fontsize=18)
ax[0].set_xlim([-1.7,1.7])
ax[0].set_ylim([-1.7,1.7])
ax[0].set_xticks([])
ax[0].set_yticks([])

ax[1].hlines(y=1.0,xmin=-0.3,xmax=0.3, color=lc[0],linestyle='-',linewidth=1)
ax[1].hlines(y=-1.0,xmin=-0.3,xmax=0.3, color=lc[0],linestyle='-',linewidth=1)
ax[1].vlines(x=1.0,ymin=-0.50,ymax=0.50, color=lc[0],linestyle='-',linewidth=2)
ax[1].vlines(x=-1.0,ymin=-0.50,ymax=0.50, color=lc[0],linestyle='-',linewidth=1)

ax[1].text(0.0, 0.82, r'$x$', fontsize=fs)
ax[1].plot(0,1, 'k.')
ax[1].text(0, 1.4, r'$\nabla\xi(x)$', fontsize=fs)
ax[1].arrow(0, 1, 0, 0.3, facecolor='black', width=0.01,head_width=0.05, head_length=0.1)
ax[1].axhline(y=1,xmin=0.02,xmax=0.98, color=lc[0],linestyle='--', linewidth=1)
ax[1].axvline(x=1.0,ymin=0.05,ymax=0.95, color=lc[0],linestyle='--',linewidth=1)
ax[1].text(0.82, -0.0, r'$y$', fontsize=fs)
ax[1].text(1.25, -0.3, r'$\nabla\xi(y)$', fontsize=fs)
ax[1].arrow(1, 0, 0.3, 0.0, facecolor='black', width=0.01,head_width=0.05, head_length=0.1)
ax[1].plot(1,0, 'ko')
ax[1].annotate(r'$T_x\Sigma$', xy=(1.45, 1.03), xytext=(1.5, 1.4),
        arrowprops=dict(facecolor='black', shrink=0.05,width=0.6,headwidth=8.0, headlength=4), fontsize=fs)
ax[1].text(0.82, 0.82, r'$v$', fontsize=fs)
ax[1].plot(1,1, 'k.')

t1 = np.linspace(0.3, 0.7, 100)
yy1 = [1.0 + math.exp(-0.22 / (tt - 0.295)**2) for tt in t1]
#ax[1].plot(t1,yy1,'k')
t = np.linspace(0.5, 2.7, 2000)
yy2 = [(1.0 + math.exp(-0.52 / (tt - 0.495)**2)) for tt in t]
f = interpolate.interp1d(yy2, t)
t2 = np.linspace(1.2, 1.0, 500)
yy2 = f(t2)
yy2[-1] = 0.5
#ax[1].plot(t2,yy2, 'k')

ttot = [*t1, *t2]
ytot = [*yy1, *yy2]
tck, u = interpolate.splprep([ttot, ytot], s=0.0)
unew = np.arange(0, 1.00, 0.001)
new_pts = interpolate.splev(unew, tck)
ax[1].plot(new_pts[0], new_pts[1], color=lc[0],linestyle='-',linewidth=1)
ax[1].plot(-new_pts[0], -new_pts[1], color=lc[0],linestyle='-',linewidth=1)
ax[1].plot(-new_pts[0], new_pts[1], color=lc[0],linestyle='-',linewidth=1)
ax[1].plot(new_pts[0], -new_pts[1], color=lc[0],linestyle='-',linewidth=1)

idx = int(275 / 600 * 1000)
ax[1].plot( new_pts[0][idx], new_pts[1][idx],'ko')
ax[1].plot( new_pts[0][idx], -new_pts[1][idx],'ko')

ax[1].set_xlim([-1.7,1.7])
ax[1].set_ylim([-1.7,1.7])
ax[1].set_xticks([])
ax[1].set_yticks([])

ax[1].spines['left'].set_color('none')
ax[1].spines['bottom'].set_color('none')
# Eliminate upper and right axes
ax[1].spines['right'].set_color('none')
ax[1].spines['top'].set_color('none')
# Show ticks in the left and lower axes only
ax[1].xaxis.set_ticks_position('bottom')
ax[1].yaxis.set_ticks_position('left')
ax[1].set_xlim([-1.8,1.8])
ax[1].set_ylim([-1.8,1.8])

y = np.linspace(-1.06, 1.06, 100)
x = [math.exp(-0.08 / math.fabs(yy)) * math.sin(3 * math.pi / math.fabs(yy)) for yy in y]
ax[2].plot(x,y, color=lc[0])
nz = 7
yz = [3.0 / k for k in range(3, nz)]
ax[2].plot(np.zeros(nz-3), yz, 'ko', linestyle='',  linewidth=1)
ax[2].plot(np.zeros(nz-3), -1.0 * np.array(yz), 'ko', linestyle='',  linewidth=1)

tvec = -math.cos(3*math.pi) * 3 * math.pi * math.exp(-0.08)
xc = x[-1] + y[-1] / tvec
r = math.sqrt((xc-x[-1])**2 + y[-1] ** 2)
theta = math.atan2(y[-1], x[-1]-xc)
t = np.linspace(0, theta, 100)
ax[2].plot(xc + r * np.cos(t), r * np.sin(t), color=lc[0])
ax[2].plot(xc + r * np.cos(t), -r * np.sin(t), color=lc[0])
ax[2].axhline(y=r,xmin=0.15,xmax=0.98, color=lc[0],linestyle='--', linewidth=1)
ax[2].text(xc+0.1, r+0.3, r'$\nabla\xi(x)$', fontsize=fs)
ax[2].annotate(r'$T_x\Sigma$', xy=(1.6, r+0.05), xytext=(1.7, r+0.3),
        arrowprops=dict(facecolor='black', shrink=0.05,width=0.6,headwidth=8.0, headlength=4), fontsize=fs)
ax[2].text(xc, r-0.15, r'$x$', fontsize=fs)
ax[2].arrow(xc, r, 0, 0.25, facecolor='black', width=0.01,head_width=0.10, head_length=0.1)
ax[2].axvline(x=0,ymin=0.05,ymax=0.95, color=lc[0],linestyle='--',linewidth=1)
ax[2].plot(0,r, 'k.')
ax[2].plot(xc,r, 'k.')
ax[2].text(0.1, r+0.04, r'$v$', fontsize=fs)
ax[2].arrow(0, 0, -0.5, 0.0, facecolor='black', width=0.01,head_width=0.05, head_length=0.1)
ax[2].text(-1.3, 0, r'$\nabla\xi(y)$', fontsize=fs)
ax[2].text(0.1, -0.05, r'$y$', fontsize=fs)

#ax[2].spines['left'].set_position('center')
ax[2].spines['left'].set_color('none')
ax[2].spines['bottom'].set_color('none')
# Eliminate upper and right axes
ax[2].spines['right'].set_color('none')
ax[2].spines['top'].set_color('none')
# Show ticks in the left and lower axes only
ax[2].xaxis.set_ticks_position('bottom')
ax[2].yaxis.set_ticks_position('left')
ax[2].set_xlim([-2.0,2.0])
ax[2].set_ylim([-1.4,1.5])
ax[2].set_xticks([])
ax[2].set_yticks([])

plt.yticks(fontsize=18)
fig.tight_layout()
out_fig_name = './illustrative.eps' 
fig.savefig(out_fig_name)

