# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 22:14:02 2020

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import copy

n = 64
map_size = 6
unit_size = map_size/n

time = 100
dt = 0.1

x = np.linspace(-map_size/2, map_size/2, n)
y = np.linspace(-map_size/2, map_size/2, n)
X,Y = np.meshgrid(x, y)

T = [[0 for i in range(n)] for j in range(n)]
for i in range(n):
    for j in range(n):
        if (i-32)**2+(j-32)**2 < 100:
            T[i][j] = 100
        if (i-32)**2+(j-63)**2 < 400:
            T[i][j] = None
        if (i-32)**2+j**2 < 400:
            T[i][j] = None
            
T = (1 - X / 2 + X**5 + Y**3) * np.exp(-X**2 - Y**2)
            
n=64
def Grad(field):
    Gradient = [[[0 for i in range(n)] for j in range(n)]for k in range(2)]
    for i in range(n):
        for j in range(n):
            if field[i][j] != None:
                G = 0
                num = 0
                try:
                    G += field[i][j]-field[i][j-1]
                    num += 1
                except:
                    pass
                    
                try:
                    G += field[i][j+1]-field[i][j]
                    num += 1
                except:
                    pass
                Gradient[0][i][j]=G/(num*(unit_size))
                
                G = 0
                num = 0
                try:
                    G += field[i][j]-field[i-1][j]
                    num += 1
                except:
                    pass
                
                try:
                    G += field[i+1][j]-field[i][j]
                    num += 1
                except:
                    pass
                Gradient[1][i][j] = G/(num*(unit_size))
            else:
                Gradient[0][i][j] = np.nan
                Gradient[1][i][j] = np.nan
                
    return Gradient

G_T = np.array(Grad(T))

level = np.linspace(T.reshape(-1,1).min(), T.reshape(-1,1).max(), 50)

fig = plt.figure(figsize=(10, 8))

#  Varying density along a streamline
ax0 = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
cp = ax0.contourf(X, Y, T, 8, levels = level, cmap = 'hot')
ax0.contour(X, Y, T, 8, colors='black',levels = level, linewidth=.5)
fig.colorbar(cp)  
ax0.streamplot(X, Y, G_T[0], G_T[1], density=[.5, 1])
ax0.set_title('Varying Density')
plt.show()