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
    areas = np.array([])
    def __init__(self,A_func,lst:list):
        # A_func is a function of x
        # ex: def A(x): return (4+x)
        self.area = A_func
        self.x0 = lst[0]
        self.x1 = lst[1]
    def plot_nozzle(self):
        x = np.linspace(self.x0,self.x1,1000*(self.x1-self.x0))
        y = self.area(x)
        plt.figure()
        plt.plot(x,y,color='blue')
        plt.plot(x,-y,color='blue')
        plt.show()

def area(x):
    y = 1+2.2*(x-1.5)**2
    return y
nozzle = nozzle(area,[-1,4])
nozzle.plot_nozzle()


