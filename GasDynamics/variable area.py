import numpy as np
from scipy.optimize import root_scalar, fsolve

"""Analyze compressible flows in nozzles accounting only for tha effect of area changes"""
class flow(object):
    def __init__(self,R=287.0,gamma=1.4):
        self.R = R
        self.gamma = gamma

    #variable_area functions
    def M2_M1(self,unknown,mach1,mach2,A2_A1,ds=0,x0=0.001,x1=0.002):
        """this function solves the equation that relates M1 TO M2 using the area ratio

            args:
            - unknown: the unknown variable in the equation, it can be set to 1, 2, or 3
                        1 : solve for mach1
                        2: solve for mach2
                        3: solve for A2_A1
            - x0 and x1 are the initial guesses used to solve the equation, they should be modified for what suits each case.
        """

        if unknown not in [1, 2, 3]:
            raise ValueError("Unknown variable must be 1, 2, or 3")
        def equation(x):
            """
            The isentropic flow equation to be solved.
            """
            if unknown == 1:
                return (
                        A2_A1 - (x / mach2) * (
                        (1 + 0.5 * (self.gamma - 1) * mach2 ** 2) / (1 + 0.5 * (self.gamma - 1) * x ** 2)
                ) ** ((self.gamma + 1) / (2 * (self.gamma - 1))) * np.exp(-ds / self.R)
                )
            elif unknown == 2:
                return (
                        A2_A1 - (mach1 / x) * (
                        (1 + 0.5 * (self.gamma - 1) * x ** 2) / (1 + 0.5 * (self.gamma - 1) * mach1 ** 2)
                ) ** ((self.gamma + 1) / (2 * (self.gamma - 1))) * np.exp(-ds / self.R)
                )
            else:
                return (
                        x - (mach1 / mach2) * (
                        (1 + 0.5 * (self.gamma - 1) * mach2 ** 2) / (1 + 0.5 * (self.gamma - 1) * mach1 ** 2)
                ) ** ((self.gamma + 1) / (2 * (self.gamma - 1))) * np.exp(-ds / self.R)
                )

        solution = root_scalar(equation, x0=x0,x1=x1)
        return solution.root

    #----------------------------------------------------------------
    # To calculate the local pressure and temperature ratios (stagnation/static)
    def mach_pressure(self,unknown, mach, pressure_ratio,x0=0.001,x1=0.002):
        '''Used to find Mach number for given stagnation pressure and gamma'''
        if unknown not in [1, 2]:
            raise ValueError("Unknown variable must be 1, 2, or 3")
        def equation(x):
            if unknown == 1:
                return (pressure_ratio - (1 + 0.5 * (self.gamma - 1) * x ** 2) ** (-self.gamma / (self.gamma - 1)))
            else:
                return (x - (1 + 0.5 * (self.gamma - 1) * mach ** 2) ** (-self.gamma / (self.gamma - 1)))

        solution = root_scalar(equation, x0=x0, x1=x1)
        return solution.root

    def mach_temp(self, unknown, mach, temp_ratio, x0=0.001,x1=0.002):
        '''Used to find Mach number for given stagnation pressure and gamma'''
        if unknown not in [1, 2]:
            raise ValueError("Unknown variable must be 1, 2, or 3")
        def equation(x):
            if unknown == 1:
                return (temp_ratio - (1 + 0.5 * (self.gamma - 1) * x ** 2))
            else:
                return (x - (1 + 0.5 * (self.gamma - 1) * mach ** 2))

        solution = root_scalar(equation, x0=x0, x1=x1)
        return solution.root

    #-----------------------------Convergant or divergant nozzle-----------------------------------
    def get_Pstar(self,pt):
        return pt * self.mach_pressure(2,mach=1,pressure_ratio=None)

    # CD NOZZLE
    def mach_ref_area(self,unknown, mach, area_ratio):
        '''Used to find Mach number for given reference area ratio and gamma'''
        if unknown not in [1, 2]:
            raise ValueError("Unknown variable must be 1, or 2 ")
        def equation(x):
            if unknown==1:
                return (area_ratio - ((1.0 / x) * ((1 + 0.5 * (self.gamma - 1) * x ** 2) /
                    ((self.gamma + 1) / 2)) ** ((self.gamma + 1) / (2 * (self.gamma - 1)))
                ))
            elif unknown==2:
                return (x - ((1.0 / mach) * ((1 + 0.5*(self.gamma - 1) * mach ** 2) /
                    ((self.gamma + 1) / 2)) ** ((self.gamma + 1) / (2 * (self.gamma - 1)))
                ))
    def get_Pcr3(self,Me,pt):
        pcr3 = pt * self.mach_pressure(unknown=2,mach=Me,pressure_ratio=None)
        return pcr3
    def get_pcr1(self,Me,pt,x0=0.001,x1=0.002):
        area_ratio = self.M2_M1(unknown=3, mach1=1, mach2=Me, A2_A1=None)
        me2 = self.M2_M1(unknown=2,mach1=1,mach2=None,A2_A1=area_ratio,x0=x0,x1=x1)
        pcr1 = pt * self.mach_pressure(unknown=2, mach=me2,pressure_ratio=None)
        return pcr1


flow1 = flow()
#Testing
print(flow1.M2_M1(unknown=2,mach1=0.5,A2_A1=2.5,mach2=None))

print(flow1.M2_M1(unknown=1,mach2=0.1759983817314642
,A2_A1=2.5,mach1=None))

print(flow1.M2_M1(unknown=3,mach2=0.1759983817314642
,A2_A1=None,mach1=0.5))

print(flow1.get_Pstar(680))

print(flow1.get_Pcr3(1.75,1000))

print(flow1.get_pcr1(1.75,1000))