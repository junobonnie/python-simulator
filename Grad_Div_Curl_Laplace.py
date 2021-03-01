# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 02:01:50 2020

@author: user
"""

import numpy as np

class BoundaryConditionError(Exception):
    """경계조건과 관련된 예외"""
    pass

class __GDCL__:
    def __init__(self, n, map_size, UP='close', DOWN='close', LEFT='close', RIGHT='close'):
        self.n = n
        self.map_size = map_size
        self.unit_size = map_size/n
        self.UP = UP
        self.DOWN = DOWN
        self.LEFT = LEFT
        self.RIGHT = RIGHT
        boundary_list = [self.UP, self.DOWN, self.LEFT, self.RIGHT]
        boundary_condition_list =['close', 'open', 'loop']
        for bl in boundary_list:
            Error = True
            for bcl in boundary_condition_list:
                if bl==bcl:
                    Error = False
                    break
            if Error:
                raise BoundaryConditionError('''경계조건 값은 'close', 'open', 'loop'만 가능''')
                
        if (self.UP=='loop' and self.DOWN!='loop') or (self.UP!='loop' and self.DOWN=='loop'):
            raise BoundaryConditionError('''UP이 'loop'면 DOWN도 'loop'여야함. 그 반대도 마찬가지임''')
            
        elif (self.LEFT=='loop' and self.RIGHT!='loop') or (self.LEFT!='loop' and self.RIGHT=='loop'):
            raise BoundaryConditionError('''LEFT가 'loop'면 RIGHT도 'loop'여야함. 그 반대도 마찬가지임''') 
                                   
    def Grad(self, field):
        Gradient = [[[0.0 for i in range(self.n)] for j in range(self.n)]for k in range(2)]
        for i in range(self.n):
            for j in range(self.n):
                if field[i][j] != None:
                    G = 0
                    try:
                        if j>0:
                            G += field[i][j]-field[i][j-1]
                    except:
                        pass
                        
                    try:
                        G += field[i][j+1]-field[i][j]
                    except:
                        pass
                    Gradient[0][i][j] = G/(2*self.unit_size)
                    
                    G = 0
                    try:
                        if i>0:
                            G += field[i][j]-field[i-1][j]
                    except:
                        pass
                    
                    try:
                        G += field[i+1][j]-field[i][j]
                    except:
                        pass
                    Gradient[1][i][j] = G/(2*self.unit_size)
                else:
                    Gradient[0][i][j] = None
                    Gradient[1][i][j] = None
                    
        return Gradient
    
    def Div(self, vector_field):
        Divergence = [[0.0 for i in range(self.n)] for j in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                if vector_field[0][i][j] != None:
                    D = 0
                    try:
                        if j>0:
                            D += vector_field[0][i][j]-vector_field[0][i][j-1]
                        else:
                            if self.LEFT == 'close':
                                D += vector_field[0][i][j]
                    except:
                        if self.LEFT == 'close':
                            D += vector_field[0][i][j]
                        pass
                        
                    try:
                        D += vector_field[0][i][j+1]-vector_field[0][i][j]
                    except:
                        if self.RIGHT == 'close':
                            D += -vector_field[0][i][j]
                        pass
                    Divergence[i][j] += D/(2*self.unit_size)
                    
                    D = 0
                    try:
                        if i>0:
                            D += vector_field[1][i][j]-vector_field[1][i-1][j]
                        else:
                            if self.DOWN == 'close':
                                D += vector_field[1][i][j]
                    except:
                        if self.DOWN == 'close':
                            D += vector_field[1][i][j]
                        pass
                    
                    try:
                        D += vector_field[1][i+1][j]-vector_field[1][i][j]
                    except:
                        if self.UP == 'close':
                            D += -vector_field[1][i][j]
                        pass
                    Divergence[i][j] += D/(2*self.unit_size)
                else:
                    Divergence[i][j] = None
                    
        return Divergence
    
    def Curl(self, vector_field):
        curl = [[0.0 for i in range(self.n)] for j in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                if vector_field[0][i][j] != None:
                    C = 0
                    try:
                        if j>0:
                            C += vector_field[1][i][j]-vector_field[1][i][j-1]
                        #else:
                        #    C += vector_field[1][i][j]
                    except:
                        #C += vector_field[1][i][j]
                        pass
                        
                    try:
                        C += vector_field[1][i][j+1]-vector_field[1][i][j]
                    except:
                        #C += -vector_field[1][i][j]
                        pass

                    try:
                        if i>0:
                            C -= vector_field[0][i][j]-vector_field[0][i-1][j]
                        #else:
                        #    C -= vector_field[0][i][j]
                    except:
                        #C -= vector_field[0][i][j]
                        pass
                    
                    try:
                        C -= vector_field[0][i+1][j]-vector_field[0][i][j]
                    except:
                        #C -= -vector_field[0][i][j]
                        pass
                    curl[i][j] += C/(2*self.unit_size)
                else:
                    curl[i][j] = None
                    
        return curl
    
    def Laplacian(self, field):
        Laplace = [[0.0 for i in range(self.n)] for j in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                L = 0
                num = 0
                try:
                    if i>0:
                        L += field[i-1][j]
                        num += 1
                except:
                    pass
                
                try:
                    L += field[i+1][j]
                    num += 1
                except:
                    pass
                
                try:
                    if j>0:
                        L += field[i][j-1]
                        num += 1
                except:
                    pass
                
                try:
                    L += field[i][j+1]
                    num += 1
                except:
                    pass
                
                try:
                    L -= num*field[i][j]
                    Laplace[i][j] = L/(self.unit_size**2)
                except:
                    Laplace[i][j] = None
                    
        return Laplace

if __name__ == "__main__": 
    import matplotlib.pyplot as plt
    from matplotlib import animation
    import copy
    
    n = 64
    map_size = 6
    unit_size = map_size/n
    GDCL = __GDCL__(n, map_size) 
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
    G_T = np.array(GDCL.Grad(T))
    
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
    
    
    
    
    D_G_T = np.array(GDCL.Div(G_T))
    
    level = np.linspace(D_G_T.reshape(-1,1).min(), D_G_T.reshape(-1,1).max(), 50)
    fig1 = plt.figure(figsize=(10, 8))
    
    #  Varying density along a streamline
    ax1 = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
    cp = ax1.contourf(X, Y, D_G_T, 8, levels = level, cmap = 'hot')
    ax1.contour(X, Y, D_G_T, 8, colors='black',levels = level, linewidth=.5)
    fig1.colorbar(cp)  
    ax1.set_title('Varying Density')
    plt.show()
    
    L_T = np.array(GDCL.Laplacian(T))
    
    level = np.linspace(L_T.reshape(-1,1).min(), L_T.reshape(-1,1).max(), 50)
    fig2 = plt.figure(figsize=(10, 8))
    
    #  Varying density along a streamline
    ax2 = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
    cp = ax2.contourf(X, Y, L_T, 8, levels = level, cmap = 'hot')
    ax2.contour(X, Y, L_T, 8, colors='black',levels = level, linewidth=.5)
    fig2.colorbar(cp)  
    ax2.set_title('Varying Density')
    plt.show()
    
    
    
    C_G_T = np.array(GDCL.Curl(G_T))
    
    level = np.linspace(C_G_T.reshape(-1,1).min(), C_G_T.reshape(-1,1).max(), 50)
    fig3 = plt.figure(figsize=(10, 8))
    
    #  Varying density along a streamline
    ax3 = plt.axes(xlim=(-map_size/2, map_size/2), ylim=(-map_size/2, map_size/2))
    cp = ax3.contourf(X, Y, C_G_T, 8, levels = level, cmap = 'hot')
    ax3.contour(X, Y, C_G_T, 8, colors='black',levels = level, linewidth=.5)
    fig3.colorbar(cp)  
    ax3.set_title('Varying Density')
    plt.show()
    
    GDCL2 = __GDCL__(n, map_size,RIGHT='loop') 