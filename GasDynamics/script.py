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
grid = 100000
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
            dA2_dx[i] = 6*c3*(x[i]-x1)

        elif x[i]>x2:
            A[i] = d1 + d2*(x[i]-x2)**2 + d3*(x[i]-x2)**3
            dA_dx[i] = 2*d2*(x[i]-x2) + 3*d3*(x[i]-x2)**2
            dA2_dx[i] = 6*d3*(x[i]-x2)
    D = np.sqrt(4*A/np.pi)
    return [A,dA_dx,dA2_dx, D]

def dm_dx(a,b,c,x,M):
    roots = np.roots([a,b,c])
    if ((x < x2) and (M<1)) or ((x>x2) and (M>1)):
        return np.max(roots)
    else:
        return np.min(roots)

def solve(f,Q,TTin,PTin,M_in,grid):
    X = np.linspace(0,1,grid+1)
    A,dA_dx,dA2_dx,D = areas(X)
    dD_dx = dA_dx/(8*np.pi*D)
    dx = X[2] - X[1]
    Ma = np.append(M_in, np.zeros(grid))
    TT = np.append(TTin, np.zeros(grid))
    PT = np.append(PTin, np.zeros(grid))
    RT = np.append(Dtin, np.zeros(grid))
    TS = np.append(Ts(Ttin, M_in), np.zeros(grid))
    PS = np.append(Ps(Ptin, M_in), np.zeros(grid))
    RS = np.append((PS[0]/(R * TS[0])), np.zeros(grid))
    SSin = np.sqrt(gamma * R * TS[0])
    SS = np.append(SSin,np.zeros(grid))
    Vin = M_in*SSin
    V = np.append(Vin, np.zeros(grid))
    dq = Q * dx / length
    a = 0
    b = 0
    c = 0
    for i in range(len(X)-1):
        Vx = V[i]
        Mx = Ma[i]
        dAx = A[i+1] - A[i]
        Dx = D[i]
        TTx = TT[i]
        XM = (1 + 0.5 * (gamma - 1) * Mx ** 2)
        XX = (Mx ** 2 - 1)
        dDx = dD_dx[i]
        dA2 = dA2_dx[i]
        if np.abs(Mx - 1) < 0.00001:
            print(f"WARNING: M = 1 at x = {X[i]}")
            break
        else:
            dM_M = (((1 + gamma * Mx ** 2) * XM / (2 * (1 - Mx ** 2))) * (dq / cp / TTx)
                    + (gamma * XM * Mx ** 2 / (2 * (1 - Mx ** 2))) * (f * dx / Dx)
                    - ((2 + (gamma - 1) * Mx ** 2) / (2 * (1 - Mx ** 2))) * (dAx / A[i]))
        dM = dM_M * Mx
        dV_V = 0.5 * (dq / cp / TTx) + (1 / XM) * dM_M
        dV = dV_V * Vx
        V[i + 1] = Vx + dV
        Ro = RS[i]
        dRo = - Ro * (dV / Vx + dAx / A[i])
        RS[i + 1] = Ro + dRo
        T = TS[i]
        dT = - T * ((gamma - 1) * (Mx ** 2) * (dV / Vx) - XM * (dq / cp / TTx))
        TS[i + 1] = T + dT
        P = PS[i]
        dP = - P * (0.5 * gamma * (Mx ** 2) * (f * dx / Dx) + gamma * (Mx ** 2) * (dV / Vx))
        PS[i + 1] = P + dP
        SS[i + 1] = np.sqrt(gamma * R * TS[i + 1])
        Ma[i + 1] = Mx + dM
        TT[i + 1] = TS[i + 1] * (1 + 0.5 * (gamma - 1) * Ma[i + 1] ** 2)
        PT[i + 1] = PS[i + 1] * ((1 + 0.5 * (gamma - 1) * Ma[i + 1] ** 2) ** (gamma / (gamma - 1)))
        RT[i + 1] = PT[i + 1] / (R * TT[i + 1])
        a = 8 / (gamma + 1)
        b = 2 * gamma / cp / TTx * dq / dx + 2 * gamma * f / Dx
        c = gamma * f * (-dDx / Dx) + 2 * (dAx) ** 2 / A[i] - 2 * dA2 / A[i] - (1 + gamma) / (cp * TTx ** 2) * dq / dx * dT / dx
        if i == 88033:
            print(dM / dx)
    return Ma,V,TT,TS,PT,PS,RT,RS

X = np.linspace(0,1,grid+1)
Ma1,V1,TT1,TS1,PT1,PS1,RT1,RS1 = solve(f,Q,Ttin,Ptin,M_in,grid)
Ma2,V2,TT2,TS2,PT2,PS2,RT2,RS2 = solve(f,0,Ttin,Ptin,M_in,grid)
Ma3,V3,TT3,TS3,PT3,PS3,RT3,RS3 = solve(0,0,Ttin,Ptin,M_in,grid)

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
plt.title(r'Mach number along the x')
plt.legend()
fig1.savefig('mach number', dpi=500)

fig2 = plt.figure()
plt.plot(X,TT1,label='area change + friction + Heat addition')
plt.plot(X,TT2,label='area change + friction')
plt.plot(X,TT3,label='area change')
plt.grid(True)
plt.title(r'Total Temperature')
plt.legend()
fig2.savefig('total temperature', dpi=500)

fig3 = plt.figure()
plt.plot(X,PT1,label='area change + friction + Heat addition')
plt.plot(X,PT2,label='area change + friction')
plt.plot(X,PT3,label='area change')
plt.grid(True)
plt.title(r'Total Pressure')
plt.legend()
fig3.savefig('pressure',dpi=500)

fig4 = plt.figure()
plt.plot(X,TS1,label='area change + friction + Heat addition')
plt.plot(X,TS2,label='area change + friction')
plt.plot(X,TS3,label='area change')
plt.grid(True)
plt.title('static Temperature')
plt.legend()
fig4.savefig('static temperature',dpi=500)

fig5 = plt.figure()
plt.plot(X,PS1,label='area change + friction + Heat addition')
plt.plot(X,PS2,label='area change + friction')
plt.plot(X,PS3,label='area change')
plt.grid(True)
plt.title('static pressure')
plt.legend()
fig5.savefig('static_pressure',dpi=500)

fig6 = plt.figure()
plt.plot(X,V1,label='area change + friction + Heat addition')
plt.plot(X,V2,label='area change + friction')
plt.plot(X,V3,label='area change')
plt.grid(True)
plt.title(r'Velocity')
plt.legend()
fig6.savefig('velocity', dpi=500)
A,dA_dx,dA2_dx,D = areas(X)
figg = plt.figure()
plt.plot(X,D,color='blue')
plt.plot(X, -1 * D, color='blue')
plt.grid(True)
plt.title("Diameter")
figg.savefig('diameter',dpi=500)
plt.show()



