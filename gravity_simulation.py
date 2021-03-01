# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 10:12:50 2020

@author: user
"""
import math as m
import random as r
import numpy as np
import matplotlib.pyplot as plt

planets = []#[{'mass':24,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]},
           #{'mass':0.024,'position':[24,0],'velocity':[0,1],'acceleration':[0,0]}]

plus_minus=[-1,+1]
max_mass=1
max_position=40
max_velocity=10.0
for i in range(10):
    planets.append({'mass':0,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})
    planets[i]['mass'] =max_mass*r.random()
    planets[i]['position'] = [r.choice(plus_minus)*max_position*r.random(),r.choice(plus_minus)*max_position*r.random()]
    planets[i]['velocity'] = [r.choice(plus_minus)*max_velocity*r.random(),r.choice(plus_minus)*max_velocity*r.random()]
planets.append({'mass':1000,'position':[0,0],'velocity':[0,0],'acceleration':[0,0]})

records = []
for planet in planets:    
    records.append({'pos_x':[planet['position'][0]],'pos_y':[planet['position'][1]]})
    
dt = 0.01
G = 1
time = 1000
def distance(pos1,pos2):
    D = m.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)
    return D

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

def recorded():
    for i in range(len(planets)):    
        records[i]['pos_x'].append(planets[i]['position'][0])
        records[i]['pos_y'].append(planets[i]['position'][1])
        
for i in range(int(time/dt)): 
    print(i/int(time/dt)*100)     
    gravity()
    velocity()
    position()
    recorded()
    
size = 40*2
plt.figure(figsize=(9, 9)) 
for record in records:
    plt.plot(record['pos_x'],record['pos_y'],'-')
    
plt.xlim([-size, size])
plt.ylim([-size, size])
plt.show()

t = np.arange(0,time+dt,dt)

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(9, 9))
ax = fig.gca(projection='3d')
for record in records:
    ax.plot(record['pos_x'],record['pos_y'],t,'-')

plt.xlim([-size, size])
plt.ylim([-size, size])
plt.show()

