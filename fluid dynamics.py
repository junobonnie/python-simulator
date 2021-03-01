# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 22:45:23 2020

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt
import Grad_Div_Curl_Laplace2 as gdcl
from matplotlib import animation
import copy

from numba import jit

n = 64
map_size = 10
unit_size = map_size/n
GDCL = gdcl.__GDCL__(n, map_size,)# LEFT='loop', RIGHT='loop')
time = 100
dt =  0.1

x = np.linspace(-map_size/2, map_size/2, n)
y = np.linspace(-map_size/2, map_size/2, n)
X,Y = np.meshgrid(x, y)

ro = np.array([[10.0 for i in range(n)] for j in range(n)], dtype=np.float64)
P = np.array([[1.0 for i in range(n)] for j in range(n)], dtype=np.float64)
v = [np.array([[0.0 for i in range(n)] for j in range(n)], dtype=np.float64), np.array([[0.0 for i in range(n)] for j in range(n)], dtype=np.float64)]

for i in range(n):
    for j in range(n):
        if (i-32)**2 < 100:
            if (j-32)**2 < 400:
                #v[0][i][j] = 0.0003
                ro[i][j] = 11
# for i in range(n):
#     for j in range(n):
#         v[0][i][j] = 0.03*(i-32)/32
#         v[1][i][j] = -0.03*(j-32)/32
# for i in range(n):
#     v[0][i][0] = 0.03*10
#     v[0][i][1] = 0.03*10
#     v[0][5][1] = 0.1
#     for j in range(n):
#         if (i-32)**2+(j-16)**2 < 25:
#             v[0][i][j] = np.nan
#             ro[i][j] = np.nan
#             P[i][j] = np.nan
eta=0.1*0
ro_list = []
v_list = []

for i in range(int(time/dt)):
    print(i/int(time/dt)*100)
    ro_list.append(copy.deepcopy(ro))
    v_list.append(copy.deepcopy(v))    
    u=v
    v[0] = u[0] + (-u[0]*np.array(GDCL.Grad(u[0])[0], dtype=np.float64)-u[1]*np.array(GDCL.Grad(u[0])[1], dtype=np.float64)-np.array(GDCL.Grad(P)[0], dtype=np.float64)/ro+eta*np.array(GDCL.Laplacian(u[0]), dtype=np.float64)/ro+eta*np.array(GDCL.Grad(np.array(GDCL.Div(u)))[0])/ro/3)*dt
    v[1] = u[1] + (-u[0]*np.array(GDCL.Grad(u[1])[0], dtype=np.float64)-u[1]*np.array(GDCL.Grad(u[1])[1], dtype=np.float64)-np.array(GDCL.Grad(P)[1], dtype=np.float64)/ro+eta*np.array(GDCL.Laplacian(u[1]), dtype=np.float64)/ro+eta*np.array(GDCL.Grad(np.array(GDCL.Div(u)))[1])/ro/3)*dt
    ro = ro - np.array(GDCL.Div([ro*u[0],ro*u[1]]), dtype=np.float64)*dt
    P =  1 + (ro - 10)/10
    # for i in range(n):
    #     v[0][i][0] = 0.03*10
    #     v[0][i][1] = 0.03*10
    #     ro[i][0] = 11
    #     for j in range(n):
    #         if ro[i][j] < 0.01:
    #             ro[i][j] = 0.01
#'''
pure_ro_list = []
for t in range(time):
    for i in range(n):
        for j in range(n):
            if not np.isnan(ro_list[t*int(1/dt)][i][j]):
                pure_ro_list.append(ro_list[t*int(1/dt)][i][j])
level = np.linspace(min(pure_ro_list), max(pure_ro_list), 50)
fig = plt.figure(figsize=(10, 8))
    
#  Varying density along a streamline
ax = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
cp = ax.contourf(X, Y, ro_list[0], 8, levels = level, cmap = 'plasma')
#ax.contour(X, Y, ro, 8, colors='black',levels = level, linewidth=.5)
cbar = fig.colorbar(cp)
#cbar.set_label('$\rho$', fontsize=18)  
#ax.streamplot(X, Y, v[0], v[1], density=[.5, 1])
ax.set_title('Density 0.0s', fontsize=20)

def animate(i):
    ax.clear()
    ax.set_title('Density %.1fs'%i, fontsize=20)
    ax.contourf(X, Y, ro_list[i*int(1/dt)], 8, levels = level, cmap = 'plasma')
    #C=ax.contour(X, Y, ro_list[i*int(1/dt)], 8, colors='black',levels = level, linewidth=.5)
    #ax.streamplot(X, Y, v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1], density=[.5, 1])
    # ax2.clabel(C, inline=1, fontsize=10) 

#call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate,
                      frames=time, interval=10,
                      blit=False)
anim.save('dynamic_images6.gif',fps=60)

plt.show()

#'''

fig2 = plt.figure(figsize=(10, 8))

ax2 = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
M_list=[]
for i in range(time):
    M = np.hypot(v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1])
    M_list.append(M)
pure_M_list = []
for t in range(time):
    for i in range(n):
        for j in range(n):
            if not np.isnan(M_list[t][i][j]):
                pure_M_list.append(M_list[t][i][j])
 
norm = plt.cm.colors.Normalize(vmin=min(pure_M_list), vmax=max(pure_M_list))
cp = ax2.quiver(X, Y, v_list[0][0], v_list[0][1], M_list[0], norm=norm, units='x', width=0.022, scale=1 / 30)
fig2.colorbar(cp)
#ax.contour(X, Y, ro, 8, colors='black',levels = level, linewidth=.5)
#ax2.streamplot(X, Y, v[0], v[1], density=[.5, 1], color='r')
ax2.set_title('Velocity 0.0s', fontsize=20)

def animate2(i):
    ax2.clear()
    ax2.set_title('Velocity %.1fs'%i, fontsize=20)
    ax2.quiver(X, Y, v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1], M_list[i], norm=norm, units='x', width=0.022, scale=1 / 30)
    #C=ax.contour(X, Y, ro_list[i*int(1/dt)], 8, colors='black',levels = level, linewidth=.5)
    #ax2.streamplot(X, Y, v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1], density=[.5, 1], color='r')
    # ax2.clabel(C, inline=1, fontsize=10) 

#call the animator.  blit=True means only re-draw the parts that have changed.
anim2 = animation.FuncAnimation(fig2, animate2,
                      frames=time, interval=10,
                      blit=False)
anim2.save('dynamic_images7.gif',fps=60)

plt.show()



# fig3 = plt.figure(figsize=(10, 8))

# ax3 = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
# cp = ax3.contourf(X, Y, ro_list[0], 8, levels = level, cmap = 'plasma')
# #ax.contour(X, Y, ro, 8, colors='black',levels = level, linewidth=.5)
# fig3.colorbar(cp)
# M = np.hypot(v_list[0][0], v_list[0][1])
# ax3.quiver(X, Y, v_list[0][0], v_list[0][1], 10*M, units='x', pivot='tip', width=0.022, scale=1 / 30)
# #ax.contour(X, Y, ro, 8, colors='black',levels = level, linewidth=.5)
# #ax2.streamplot(X, Y, v[0], v[1], density=[.5, 1], color='r')
# ax3.set_title('Density & Velocity 0.0s', fontsize=20)

# def animate3(i):
#     ax3.clear()
#     ax3.set_title('Density & Velocity %.1fs'%i, fontsize=20)
#     ax3.contourf(X, Y, ro_list[i*int(1/dt)], 8, levels = level, cmap = 'plasma')
#     M = np.hypot(v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1])
#     ax3.quiver(X, Y, v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1], 10*M, units='x', pivot='tip', width=0.022, scale=1 / 30)
#     #C=ax.contour(X, Y, ro_list[i*int(1/dt)], 8, colors='black',levels = level, linewidth=.5)
#     #ax2.streamplot(X, Y, v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1], density=[.5, 1], color='r')
#     # ax2.clabel(C, inline=1, fontsize=10) 

# #call the animator.  blit=True means only re-draw the parts that have changed.
# anim3 = animation.FuncAnimation(fig3, animate3,
#                       frames=time, interval=10,
#                       blit=False)
# anim3.save('dynamic_images8.gif',fps=60)

# plt.show()
'''
fig3 = plt.figure(figsize=(12, 8))

ax3 = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
cp = ax3.contourf(X, Y, ro_list[0], 8, levels = level, cmap = 'plasma')
#ax.contour(X, Y, ro, 8, colors='black',levels = level, linewidth=.5)
fig3.colorbar(cp)
#ax2.quiver(X, Y, v_list[0][0], v_list[0][1], 10*M, units='x', pivot='tip', width=0.022, scale=1 / 30)
#ax.contour(X, Y, ro, 8, colors='black',levels = level, linewidth=.5)
cp2 = ax3.streamplot(X, Y, v_list[0][0], v_list[0][1], color=M_list[0], norm=norm, density=[1, 2], cmap='jet')
#fig3.colorbar(cp, cmap='jet')
ax3.set_title('Density & Velocity 0.0s', fontsize=20)
ax3.set_xlim(-map_size/2, map_size/2)
ax3.set_ylim(-map_size/2, map_size/2)
def animate3(i):   
    ax3.clear()   
    ax3.set_title('Density & Velocity %.1fs'%i, fontsize=20)
    ax3.contourf(X, Y, ro_list[i*int(1/dt)], 8, levels = level, cmap = 'plasma')
    #ax2.quiver(X, Y, v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1], 10*M, units='x', pivot='tip', width=0.022, scale=1 / 30)
    #C=ax.contour(X, Y, ro_list[i*int(1/dt)], 8, colors='black',levels = level, linewidth=.5)
    ax3.streamplot(X, Y, v_list[i*int(1/dt)][0], v_list[i*int(1/dt)][1], color=M_list[i], norm=norm, density=[1, 2], cmap='jet')
    # ax2.clabel(C, inline=1, fontsize=10) 
    ax3.set_xlim(-map_size/2, map_size/2)
    ax3.set_ylim(-map_size/2, map_size/2)
#call the animator.  blit=True means only re-draw the parts that have changed.
anim3 = animation.FuncAnimation(fig3, animate3,
                      frames=time, interval=10,
                      blit=False)
anim3.save('dynamic_images8.gif',fps=60)

plt.show()
'''


fig4 = plt.figure(figsize=(10, 8))

ax4 = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
C_v_list=[]
for i in range(time):
    C_v = GDCL.Curl(v_list[i*int(1/dt)])
    C_v_list.append(C_v)
C_v_list = np.array(C_v_list)

pure_C_v_list = []
for t in range(time):
    for i in range(n):
        for j in range(n):
            if not np.isnan(C_v_list[t][i][j]):
                pure_C_v_list.append(C_v_list[t][i][j])

if abs(min(pure_C_v_list)) < abs(max(pure_C_v_list)):
    C_v_max = abs(max(pure_C_v_list))
else:
    C_v_max = abs(min(pure_C_v_list))
level2 = np.linspace(-C_v_max, C_v_max, 50)
cp = ax4.contourf(X, Y, C_v_list[0], 8, levels = level2, cmap = 'seismic')
fig4.colorbar(cp)
#ax2.quiver(X, Y, v_list[0][0], v_list[0][1], 10*M, units='x', pivot='tip', width=0.022, scale=1 / 30)
#ax.contour(X, Y, ro, 8, colors='black',levels = level, linewidth=.5)
ax4.set_title('Vorticity 0.0s', fontsize=20)
def animate4(i):   
    ax4.clear()   
    ax4.set_title('Vorticity %.1fs'%i, fontsize=20)
    ax4.contourf(X, Y, C_v_list[i], 8, levels = level2, cmap = 'seismic')
    
#call the animator.  blit=True means only re-draw the parts that have changed.
anim4 = animation.FuncAnimation(fig4, animate4,
                      frames=time, interval=10,
                      blit=False)
anim4.save('dynamic_images9.gif',fps=60)

plt.show()