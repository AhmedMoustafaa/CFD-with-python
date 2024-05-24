import numpy as np
import matplotlib.pyplot as plt
"""
analyze compressible flows with multiple effect: 
    - area change
    - friction
    - heat transfer
    
    I will try to the case of shocks
"""

class nozzle(object):

    def __init__(self,A_func,dA_func,D_func,dD_func,mesh,Ttn,Ptn,M_in,Q=0,f=0,gamma=1.4,cp=1.004):
        # A_func is a function of x
        # ex: def A(x): return (4+x)
        self.area = A_func
        self.dA_dx = dA_func
        self.D = D_func
        self.dD_dx = dD_func
        self.mesh = mesh
        self.Ttn = Ttn
        self.Ptn = Ptn
        self.M_in = M_in
        self.Q = Q
        self.f = f
        self.gamma = gamma
        self.cp = cp

    

    def plot_nozzle(self):
        x = np.linspace(0,1,self.mesh)
        y = np.array([])
        for i in range(len(x)):
            y = np.append(y, D(x[i]))
        plt.figure()
        plt.plot(x,y,color='blue')
        plt.plot(x,-y,color='blue')
        plt.show()



def areas(x):
    A=0
    dA_dx = 0
    A1 = 0.5612
    A2 = 0.1403
    A3 = 0.5612
    x1 = 0
    x2 = 0.2113
    x3 = 1
    c2 = 3*(A2-A1)/(x2-x1)**2
    c3 = -2*(A2-A1)/(x2-x1)**3
    d1 = A2
    d2 = 3*(A3-A2)/(x3-x2)**2
    d3 = -2*(A3-A2)/(x3-x2)**3
    if x<x2:
        A = A1 + c2*(x-x1)**2 + c3*(x-x1)**3
        dA_dx = 2*c2*(x-x1) + 3*c3*(x-x1)**2

    elif x>x2:
        A = d1 + d2*(x-x2)**2 + d3*(x-x2)**3
        dA_dx = 2*d2*(x-x2) + 3*d3*(x-x2)**2
    elif x==x2:
        A=A2
        dA_dx1 = 2 * c2 * (x - x1) + 3 * c3 * (x - x1) ** 2
        dA_dx2 = 2 * d2 * (x - x2) + 3 * d3 * (x - x2) ** 2
        dA_dx = (dA_dx1 + dA_dx2)/2.0
    D = A/np.pi * 4
    return [A,dA_dx, D]
def area(x):
    return areas(x)[0]
def da_dx(x):
    return areas(x)[1]
def D(x):
    return areas(x)[2]
def dD_dx(x):
    return areas(x)[1] * areas(x)[2] * 2 / np.pi

nozzle = nozzle(area,da_dx,D,dD_dx,100)
nozzle.plot_nozzle()

