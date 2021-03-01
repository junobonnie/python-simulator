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
planets = []#[{'mass':24,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]},
           #{'mass':0.024,'position':[24,0],'velocity':[0,1],'acceleration':[0,0]}]

plus_minus=[-1,+1]
max_mass=50#1
max_position=40
max_velocity=0#10.0
max_radius=20
for i in range(10):
    random_angle_1=2*m.pi*r.random()
    random_angle_2=2*m.pi*r.random()
    planets.append({'mass':0,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})
    planets[i]['mass'] =max_mass*r.random()
    planets[i]['position'] = [max_position*m.cos(random_angle_1)*r.random(),max_position*m.sin(random_angle_1)*r.random()]
    planets[i]['velocity'] = [max_velocity*m.cos(random_angle_2)*r.random(),max_velocity*m.sin(random_angle_2)*r.random()]
    planets[i]['color'] = (0.9*r.random(),0.9*r.random(),0.9*r.random())
    planets[i]['radius'] = max_radius*r.random()
#planets.append({'mass':1000,'position':[0,0],'velocity':[0,0],'acceleration':[0,0],'color':(1,0,0)})
############################################################################### 
# plus_minus=[-1,+1]
# max_mass=0
# max_position=40
# max_velocity=10.0
# for i in range(100):
#     planets.append({'mass':0,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})
#     planets[i]['mass'] =max_mass*r.random()
#     planets[i]['position'] = [-90,r.choice(plus_minus)*max_position*r.random()]
#     planets[i]['velocity'] = [7,0]
# planets.append({'mass':1000,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})
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
    records.append({'mass':[planet['mass']],'pos_x':[planet['position'][0]],'pos_y':[planet['position'][1]],'color':planet['color']})
###############################################################################    
dt = 0.01
G = 1
time = 100
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
        if planet['mass']==0:
            pass
        else:
            for i in range(2):
                planet['position'][i] += planet['velocity'][i]*dt
            

def recorded():
    for i in range(len(planets)):    
        records[i]['mass'].append(planets[i]['mass'])
        records[i]['pos_x'].append(planets[i]['position'][0])
        records[i]['pos_y'].append(planets[i]['position'][1])
        
def collision():
    for planet_1 in planets:
        for planet_2 in planets:
            if planet_1 == planet_2:
                pass
            else:
                if distance(planet_1['position'],planet_2['position']) < (planet_1['radius']+planet_2['radius']):
                    if planet_1['mass'] > planet_2['mass']:
                        planet_1['mass'] += planet_2['mass']
                        planet_2['mass'] = 0
                        
                    else:
                        planet_2['mass'] += planet_1['mass']
                        planet_1['mass'] = 0   
###############################################################################         
for i in range(int(time/dt)): 
    print(i/int(time/dt)*100)     
    gravity()
    velocity()
    position()
    recorded()
    collision()
###############################################################################   
size = 40*2
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

line=[]
point=[]
for i in range(len(planets)):
    line.append([])
    point.append([])
    
# initialization function: plot the background of each frame


# animation function.  This is called sequentially
def animate(i):
    for j in range(len(line)):
        line[i], = ax.plot(records[j]['pos_x'][:i*int(0.01*time/dt)], records[j]['pos_y'][:i*int(0.01*time/dt)], lw=2,color=(1.1*planets[j]['color'][0],1.1*planets[j]['color'][1],1.1*planets[j]['color'][2]))
        point[i], = ax.plot(records[j]['pos_x'][i*int(0.01*time/dt)], records[j]['pos_y'][i*int(0.01*time/dt)],markersize=planets[j]['radius'], marker= 'o',color=planets[j]['color'])
    
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate,
                     frames=100, interval=int(time),
                     blit=False)
anim.save('dynamic_images.gif')

plt.show()
############################################################################### 