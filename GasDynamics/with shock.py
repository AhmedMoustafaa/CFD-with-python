import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# boundry conditions
Ttin = 840.0
Ptin = 3721700.0
M_in = 0.137
Q = 400000.0
f = 0.02
gamma = 1.4
R = 287.0
cp = 1004.5
grid = 10000
length = 1.0
Dtin = Ptin/(R*Ttin)
x2 = 0.2113
def Ts(Tt,M):
    return Tt/((1 + 0.5*(gamma-1)*M**2))
def Ps(Pt,M):
    return Pt/(((1 + 0.5*(gamma-1)*M**2))**(gamma/(gamma-1)))
def areas(x):
    A=np.zeros(len(x))
    dA_dx = np.zeros(len(x))
    dA2_dx = np.zeros(len(x))
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
    for i in range(len(x)):
        if x[i]<=x2:
            A[i] = A1 + c2*(x[i]-x1)**2 + c3*(x[i]-x1)**3
            dA_dx[i] = 2*c2*(x[i]-x1) + 3*c3*(x[i]-x1)**2
            dA2_dx[i] = 2*c2 + 6*c3*(x[i]-x1)

        elif x[i]>x2:
            A[i] = d1 + d2*(x[i]-x2)**2 + d3*(x[i]-x2)**3
            dA_dx[i] = 2*d2*(x[i]-x2) + 3*d3*(x[i]-x2)**2
            dA2_dx[i] = 2*d2 + 6*d3*(x[i]-x2)
    D = np.sqrt(4*A/np.pi)
    dD_dx = 2*dA_dx / D / np.pi
    return [A,dA_dx,dA2_dx, D, dD_dx]

def dm_dx(a,b,c,x,M):
    roots = np.roots([a,b,c])
    if ((x < x2) and (M<1)) or ((x>x2) and (M>1)):
        return np.max(roots)
    else:
        return np.min(roots)

def solve(f,Q,TTin,M_in,grid):
    X = np.linspace(0,1,grid+1)
    A,dA_dx,dA2_dx,D, dD_dx = areas(X)
    dD_dx = dA_dx/(8*np.pi*D)
    dx = X[2] - X[1]
    Ma = np.append(M_in, np.zeros(grid))
    TT = np.append(TTin, np.zeros(grid))
    TS = np.append(Ts(Ttin, M_in), np.zeros(grid))
    dq = Q * dx / length
    dq_dx = dq / dx
    a = 0
    b = 0
    c = 0
    for i in range(len(X)-1):
        Mx = Ma[i]
        dAx = dA_dx[i]
        Dx = D[i]
        TTx = TT[i]
        XM = (1 + 0.5 * (gamma - 1) * Mx ** 2)
        dDx = dD_dx[i]
        dA2 = dA2_dx[i]
        a = 2 + 4/(gamma + 1)
        d = 4 * gamma ** 2 + 2 * gamma - 2 + 0.25 * (gamma**2 + 2 * gamma + 1)
        b = gamma * f / Dx + d * dq_dx / cp / TTx
        c = (dAx**2 - dA2) / A[i] - gamma * f / 2.0 / Dx * dDx
        if np.abs(Mx - 1) < 0.00001:
            dM_M = dm_dx(a,b,c,X[i],Mx)
        else:
            dM_M = (((1 + gamma * Mx ** 2) * XM / (2 * (1 - Mx ** 2))) * (dq / cp / TTx)
                    + (gamma * XM * Mx ** 2 / (2 * (1 - Mx ** 2))) * (f * dx / Dx)
                    - ((2 + (gamma - 1) * Mx ** 2) / (2 * (1 - Mx ** 2))) * (dAx / A[i]))
        T = TS[i]
        dT = T*((dq/cp/TTx) - ((gamma - 1)*Mx**2 * dM_M * XM))
        TS[i+1] = T + dT
        Ma[i + 1] = Mx + dM_M * Mx
    return Ma, TT

X = np.linspace(0,1,grid+1)
Ma1, TT1 = solve(f,Q,Ttin,M_in,grid)
Ma2, TT2 = solve(f,0,Ttin,M_in,grid)
Ma3, TT3 = solve(0,0,Ttin,M_in,grid)

print(Ma1[88033])

"""data = {'Mach1': Ma1, 'V1': V1, 'Ts1': TS1, 'TT1': TT1, 'PT1': PT1, 'PS1': PS1}
df = pd.DataFrame(data)
df.to_csv('arrays1.csv', index=False)  # Saves without the index row

data = {'Mach2': Ma2, 'V2': V2, 'Ts2': TS2, 'TT2': TT2, 'PT2': PT2, 'PS2': PS2}
df = pd.DataFrame(data)
df.to_csv('arrays2.csv', index=False)  # Saves without the index row


data = {'Mach3': Ma3, 'V3': V3, 'Ts3': TS3, 'TT3': TT3, 'PT3': PT3, 'PS3': PS3}
df = pd.DataFrame(data)
df.to_csv('arrays3.csv', index=False)  # Saves without the index row"""

fig1 = plt.figure()
plt.plot(X,Ma1,label='area change + friction + Heat addition')
plt.plot(X,Ma2,label='area change + friction')
plt.plot(X,Ma3,label='area change')
plt.grid(True)
plt.title(r'Mach number')
plt.legend()
fig1.savefig('mach number')

fig2 = plt.figure()
plt.plot(X,TT1,label='area change + friction + Heat addition')
plt.plot(X,TT2,label='area change + friction')
plt.plot(X,TT3,label='area change')
plt.grid(True)
plt.title(r'Temperature')
plt.legend()
fig2.savefig('total temperature')






