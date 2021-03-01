# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 23:15:54 2020

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import copy
import Grad_Div_Curl_Laplace as gdcl

n = 64
map_size = 10
unit_size = map_size/n

GDCL = gdcl.__GDCL__(n, map_size)

time = 100
dt = 0.1

x = np.linspace(-map_size/2, map_size/2, n)
y = np.linspace(-map_size/2, map_size/2, n)
X,Y = np.meshgrid(x, y)

T = [[0 for i in range(n)] for j in range(n)]
for i in range(n):
    for j in range(n):
        if (i-32)**2+(j-32)**2 < 100:
            T[i][j] = 0
        if (i-32)**2+(j-63)**2 < 400:
            T[i][j] = None
        if (i-32)**2+j**2 < 400:
            T[i][j] = None

e = [[0 for i in range(n)] for j in range(n)]
for i in range(n):
    for j in range(n):
        if (i-32)**2+(j-32)**2 < 100:
            e[i][j] = 1
            
Alpha = [[0.05 for i in range(n)] for j in range(n)]
Kapha = [[1 for i in range(n)] for j in range(n)]


all_T_list = []
for t in range(int(time/dt)):   
    print(t/int(time/dt)*100)
    all_T_list.append(copy.deepcopy(T))
    L_T = GDCL.Laplacian(T)
    U = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            try:
                T[i][j] += Alpha[i][j]*L_T[i][j]*dt + Kapha[i][j]*e[i][j]*dt
            except:
                T[i][j] = None
    
    
U = []
T_list = []
for t in range(time):
    for i in range(n):
        for j in range(n):
            if all_T_list[t*int(1/dt)][i][j] != None:
                U.append(all_T_list[t*int(1/dt)][i][j])
            else:
                all_T_list[t*int(1/dt)][i][j] = np.nan
    T_list.append(all_T_list[t*int(1/dt)])
    
    
level = np.linspace(min(U), max(U), 50)


fig = plt.figure(figsize=(10, 8))
ax = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
ax.set_title('Heat Simulation 0.0s', fontsize=20)
cp = ax.contourf(X, Y, all_T_list[0], 8, levels = level, cmap = 'hot')
ax.contour(X, Y, all_T_list[0], 8, colors='black',levels = level, linewidth=.5)
Alpha = np.array(Alpha)
G_T = -5*Alpha*np.array(GDCL.Grad(all_T_list[0]))
M = np.hypot(G_T[0], G_T[1])
ax.quiver(X, Y, G_T[0], G_T[1], M, units='x', pivot='tip', width=0.022, scale=1 / 0.15)# cmap='plasma')
fig.colorbar(cp)   

def animate(i):
    ax.clear()
    ax.set_title('Heat Simulation %.1fs'%i, fontsize=20)
    ax.contourf(X, Y, all_T_list[i*int(1/dt)], 8, levels = level, cmap = 'hot')
    C=ax.contour(X, Y, all_T_list[i*int(1/dt)], 8, colors='black',levels = level, linewidth=.5)
    G_T = -5*Alpha*np.array(GDCL.Grad(all_T_list[i*int(1/dt)]))
    M = np.hypot(G_T[0], G_T[1])
    ax.quiver(X, Y, G_T[0], G_T[1], M, units='x', pivot='tip', width=0.022, scale=1 / 0.15)# cmap='plasma')
    # ax2.clabel(C, inline=1, fontsize=10) 

#call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate,
                      frames=time, interval=10,
                      blit=False)
anim.save('dynamic_images3.gif',fps=60)

plt.show()

T_list = np.array(T_list)

from mpl_toolkits.mplot3d import Axes3D

fig2 = plt.figure(figsize=(10, 10))
ax2 = Axes3D(fig2)
ax2.set_title('Heat Simulation 0.0s', fontsize=30)

ax2.plot_surface(X, Y, T_list[0], rstride=1, cstride=1, vmin=min(U), vmax=max(U), cmap='hot')
ax2.contourf(X, Y, T_list[0], zdir='z', offset=min(U), levels = level, cmap='hot')
ax2.set_zlim(-min(U), max(U))

def animate2(i):
    ax2.clear()
    ax2.set_title('Heat Simulation %.1fs'%i, fontsize=30)
    ax2.plot_surface(X, Y, T_list[i], rstride=1, cstride=1, vmin=min(U), vmax=max(U), cmap='hot')
    ax2.contourf(X, Y, T_list[i], zdir='z', offset=min(U), levels = level, cmap = 'hot')
    ax2.set_zlim(-min(U), max(U))
    
#call the animator.  blit=True means only re-draw the parts that have changed.
anim2 = animation.FuncAnimation(fig2, animate2,
                      frames=time, interval=10,
                      blit=False)
anim2.save('dynamic_images4.gif',fps=60)

plt.show()

