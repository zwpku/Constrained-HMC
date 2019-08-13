#!/usr/bin/env python3
 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
 
n = 100
 
phi = np.linspace(0, 2 * np.pi, endpoint=True, num=n) 
theta = np.linspace(0, 2 * np.pi, endpoint=True, num=n) 
phi, theta = np.meshgrid(phi, theta) 
 
R=2.5 
r=1 
 
# coordinate transformation from toroidal to cartesian (the standart one)
x = np.cos(theta) * (R + r * np.cos(phi))
y = np.sin(theta) * (R + r * np.cos(phi))
z = r * np.sin(phi)
z1 = (x**2 - (R + r)**2)**2
z2 = z**2
N1 = z1 / z1.max()
N2 = z2 / z2.max()
 
fig = plt.figure(figsize=plt.figaspect(0.5), dpi=700)
#fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection='3d') 
ax2 = fig.add_subplot(1, 2, 2, projection='3d') 
plt.subplots_adjust(wspace=0.0)

ax1.plot_surface(x, y, z, facecolors=plt.cm.jet(N1), rstride=1, cstride=1, linewidth=2, shade=True)
ax2.plot_surface(x, y, z, facecolors=plt.cm.jet(N2), rstride=1, cstride=1)
 
lim = (max(abs(max(np.max(x), np.max(y), np.max(z))), abs(min(np.min(x), np.min(y), np.min(z))))) # find the absolute maximum of the function

ax1.set_xlim(-lim, lim)
ax1.set_ylim(-lim, lim)
ax1.set_zlim(-lim, lim)
ax1.view_init(elev=30, azim=25)
ax1.set_xlabel('$x$', fontsize=20, rotation=50)
ax1.set_ylabel('$y$', fontsize=20, rotation=0)
ax1.set_zlabel('$z$', fontsize=20, rotation=0)
ax1.xaxis.labelpad=-5
ax1.yaxis.labelpad=-5
ax1.zaxis.labelpad=-5
ax1.tick_params(axis='x', pad=-3) 
ax1.tick_params(axis='y', pad=-5) 
ax1.tick_params(axis='z', pad=-1) 

ax2.set_xlim(-lim, lim)
ax2.set_ylim(-lim, lim)
ax2.set_zlim(-lim, lim)
ax2.view_init(elev=30, azim=25)
ax2.set_xlabel('$x$', fontsize=20, rotation=50)
ax2.set_ylabel('$y$', fontsize=20, rotation=0)
ax2.set_zlabel('$z$', fontsize=20, rotation=0)
ax2.tick_params(axis='x', pad=-3) 
ax2.tick_params(axis='y', pad=-5) 
ax2.tick_params(axis='z', pad=-1) 
ax2.xaxis.labelpad=-5
ax2.yaxis.labelpad=-5
ax2.zaxis.labelpad=-5

# deactivae coordinate axis rendering
#plt.axis('off')
 
plt.savefig("torus.eps", bbox_inches='tight')
