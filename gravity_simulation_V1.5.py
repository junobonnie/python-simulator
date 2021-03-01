# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 15:12:15 2020

@author: user
"""

import math as m
import random as r
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

############################################################################### 
planets = []
planets = [{'mass':1000,'position':[0,0],'velocity':[0,-0.5],'acceleration':[0,0],'color':(1,0,0)},
            {'mass':100,'position':[40,0],'velocity':[0,5],'acceleration':[0,0],'color':(0.5,0,0.5)},
            {'mass':0.01,'position':[42.5,0],'velocity':[0,10],'acceleration':[0,0],'color':(0,0.1,0.8)}]

# planets = [{'mass':100,'position':[40,0],'velocity':[0,0],'acceleration':[0,0],'color':(1,0,0)},
#            {'mass':100,'position':[-40,0],'velocity':[0,-1],'acceleration':[0,0],'color':(0.5,0,0.5)},
#            {'mass':100,'position':[0,40],'velocity':[-1,0],'acceleration':[0,0],'color':(0,0.1,0.8)},
#            {'mass':100,'position':[0,-40],'velocity':[1,0],'acceleration':[0,0],'color':(0.1,0.1,0.1)}]

# planets = [{'mass':640,'position':[-0.97000436*40,0.24308753*40],'velocity':[0.4662036850*4,0.4323657300*4],'acceleration':[0,0],'color':(1,0,0)},
#             {'mass':640,'position':[0.97000436*40,-0.24308753*40],'velocity':[0.4662036850*4,0.4323657300*4],'acceleration':[0,0],'color':(0.5,0,0.5)},
#             {'mass':640,'position':[0,0],'velocity':[-0.93240737*4,-0.86473146*4],'acceleration':[0,0],'color':(0,0.1,0.8)}]
############################################################################### 
# plus_minus=[-1,+1]
# max_mass=1
# max_position=80
# max_velocity=5*m.sqrt(40)

# for i in range(100):
#     random_angle_1 = 2*m.pi*r.random()
#     random_angle_2 = 2*m.pi*r.random()
#     R = max_position*r.random()
#     v = max_velocity/m.sqrt(R)#*r.random()
#     planets.append({'mass':0,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})
#     planets[i]['mass'] = max_mass*r.random()
#     planets[i]['position'] = [R*m.cos(random_angle_1), R*m.sin(random_angle_1)]
#     planets[i]['velocity'] = [-v*m.sin(random_angle_1), v*m.cos(random_angle_1)]
#     planets[i]['color'] = (0.9*r.random(), 0.9*r.random(), 0.9*r.random())
# planets.append({'mass':1000, 'position':[0,0], 'velocity':[0,0], 'acceleration':[0,0], 'color':(1,0,0)})
############################################################################### 
# plus_minus=[-1,+1]
# max_mass=0.1
# max_position=40
# max_velocity=10.0
# for i in range(100):
#     planets.append({'mass':0,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})
#     planets[i]['mass'] =max_mass
#     planets[i]['position'] = [-90,r.choice(plus_minus)*max_position*r.random()]
#     planets[i]['velocity'] = [7,0]
#     planets[i]['color'] = (0.9*r.random(),0.9*r.random(),0.9*r.random())
# planets.append({'mass':1000,'position':[0,0],'velocity':[0,0],'acceleration':[0,0],'color':(1,0,0)})
############################################################################### 
# plus_minus=[-1,+1]
# max_mass=0.01
# max_position=10
# max_velocity=10.0
# for i in range(100):
#     planets.append({'mass':0,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})
#     planets[i]['mass'] =max_mass*r.random()
#     pm=r.choice(plus_minus)
#     planets[i]['position'] = [pm*90,r.choice(plus_minus)*max_position*r.random()]
#     planets[i]['velocity'] = [-pm*5,0]
#planets.append({'mass':1000,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})
############################################################################### 
records = []
for planet in planets:    
    records.append({'pos_x':[planet['position'][0]],'pos_y':[planet['position'][1]],'color':planet['color']})
###############################################################################    
dt = 0.1
G = 1
time = 100
size = 40*2
def distance(pos1,pos2):
    D = m.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)
    return D
############################################################################### 
def gravity():
    for planet_1 in planets:
        a = [0,0]
        for planet_2 in planets:
            if planet_1 == planet_2:
                pass
            else:
                for i in range(2):
                    a[i] += -G*planet_2['mass']/(distance(planet_1['position'],planet_2['position'])**3)*(planet_1['position'][i]-planet_2['position'][i])
                planet_1['acceleration'] = a
                
def velocity():
    for planet in planets:
        for i in range(2):
            planet['velocity'][i] += planet['acceleration'][i]*dt

def position():
    for planet in planets:
        for i in range(2):
            planet['position'][i] += planet['velocity'][i]*dt
Xmesh, Ymesh = np.meshgrid(np.linspace(-size, size, 100), np.linspace(-size, size, 100))
potential_list=[]
def potential():
    U = 0
    for planet in planets:
        U += -G*planet['mass']/np.sqrt((planet['position'][0]-Xmesh)**2+(planet['position'][1]-Ymesh)**2)
    potential_list.append(U)
def recorded():
    for i in range(len(planets)):    
        records[i]['pos_x'].append(planets[i]['position'][0])
        records[i]['pos_y'].append(planets[i]['position'][1])
###############################################################################         
for i in range(int(time/dt)): 
    print(i/int(time/dt)*100)     
    gravity()
    velocity()
    position()
    recorded()
    potential()
###############################################################################   

plt.figure(figsize=(9, 9)) 
for record in records:
    plt.plot(record['pos_x'],record['pos_y'],'-',color=record['color'])
    
plt.xlim([-size, size])
plt.ylim([-size, size])
plt.show()
############################################################################### 
t = np.arange(0,time+dt,dt)

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(9, 9))
ax = fig.gca(projection='3d')
for record in records:
    ax.plot(record['pos_x'],record['pos_y'],t,'-',color=record['color'])

plt.xlim([-size, size])
plt.ylim([-size, size])
plt.show()
############################################################################### 
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(figsize=(9, 9))
ax = plt.axes(xlim=(-size, size), ylim=(-size, size))
ax.set_title('Gravity Simulation 0.0fs', fontsize=20)
line=[]
point=[]
for i in range(len(planets)):
    line.append([])
    point.append([])
    
    line[i], = ax.plot([], [], lw=2,color=planets[i]['color'],alpha=0.7)
    point[i], = ax.plot([], [], 'o',color=planets[i]['color'])
# initialization function: plot the background of each frame
def init():
    for j in range(len(line)):
        line[j].set_data([], [])
    return line

# animation function.  This is called sequentially
def animate(i):
    ax.set_title('Gravity Simulation %.1fs'%i, fontsize=20)
    for j in range(len(line)):
        line[j].set_data(records[j]['pos_x'][:i*int(1/dt)], records[j]['pos_y'][:i*int(1/dt)])
        point[j].set_data(records[j]['pos_x'][i*int(1/dt)], records[j]['pos_y'][i*int(1/dt)])
    return line

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate,
                      frames=time, interval=10,
                      blit=True)
anim.save('dynamic_images.gif',fps=60)

plt.show()
############################################################################### 


############################################################################### 
# First set up the figure, the axis, and the plot element we want to animate

fig2 = plt.figure(figsize=(11, 9))
ax2 = plt.axes(xlim=(-size, size), ylim=(-size, size))
ax2.set_title('Gravity Simulation 0.0fs', fontsize=20)

all_potential_list = np.array([])
for i in range(time):
    all_potential_list = np.concatenate((all_potential_list,potential_list[i*int(1/dt)]), axis=None)

level = np.linspace(all_potential_list.reshape(-1, 1).min()/50, all_potential_list.reshape(-1, 1).max(), 50)
cp=ax2.contourf(Xmesh, Ymesh, potential_list[0],levels = level)#, cmap = 'plasma')
fig2.colorbar(cp)   
# initialization function: plot the background of each frame


# animation function.  This is called sequentially
def animate2(i):
    ax2.clear()
    ax2.set_title('Gravity Simulation %.1fs'%i, fontsize=20)
    ax2.contourf(Xmesh, Ymesh, potential_list[i*int(1/dt)],levels = level)#, cmap = 'plasma')
    #C=ax2.contour(Xmesh, Ymesh, potential_list[i*int(1/dt)], 20, colors='black',levels = level, linewidth=.5)
    #ax2.clabel(C, inline=1, fontsize=10)
  
    

# call the animator.  blit=True means only re-draw the parts that have changed.
anim2 = animation.FuncAnimation(fig2, animate2,
                      frames=time, interval=10,
                      blit=False)
anim2.save('dynamic_images2.gif',fps=60)

plt.show()