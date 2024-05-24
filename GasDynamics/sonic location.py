import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# boundry conditions
Ttin = 840.0
Ptin = 3721700.0
M_in = 0.2
Q = 400000
f = 0.01
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
            dA2_dx[i] = 2*c2 + 6*c3*(x[i]-x1)

        elif x[i]>x2:
            A[i] = d1 + d2*(x[i]-x2)**2 + d3*(x[i]-x2)**3
            dA_dx[i] = 2*d2*(x[i]-x2) + 3*d3*(x[i]-x2)**2
            dA2_dx[i] = 2*d2 + 6*d3*(x[i]-x2)
    D = np.sqrt(4*A/np.pi)
    dD_dx = 2*dA_dx / D / np.pi
    return [A,dA_dx,dA2_dx, D, dD_dx]


X = np.linspace(0,1,grid+1)
A, dA_dx, dA2_dx, D, dD_dx = areas(X)
min = 100
ii = 0
location = 8
dx = X[2] - X[1]
dq = Q * dx / length
dq_dx = dq / dx
for i in range(len(X)-1):
    Tx = Ttin + 398.284*X[i]
    Dx = D[i]
    dDx = dD_dx[i]
    term1 = ((gamma + 1) * gamma * f - (gamma + 1) * dDx) / Dx + 10.08 * dq_dx / 4 / cp / Tx
    if abs(term1) < 1e-4 and abs(term1) < min:
        min = abs(term1)
        location = X[i]
        ii = i
    else:
        continue
def dm_dx(a,b,c):
    roots = np.roots([a,b,c])
    return roots

D_sp = D[ii]
A_sp = A[ii]
dA_dx_sp = dA_dx[ii]
dD_dx_sp = dD_dx[ii]
dA2_dx_sp = dA2_dx[ii]
Tsp = Ttin + 398.284*X[ii]
a = 2 + 4 / (gamma + 1)
d = 4 * gamma ** 2 + 2 * gamma - 2 + 0.25 * (gamma ** 2 + 2 * gamma + 1)
b = gamma * f / D_sp + d * dq_dx / cp / Tsp
c = (dA_dx_sp ** 2 - dA2_dx_sp) / A[ii] - gamma * f / 2.0 / D_sp * dD_dx_sp
dm_dx_sp = np.max(dm_dx(a,b,c))
dm_dx_sp2 = np.min(dm_dx(a,b,c))
print(dm_dx_sp)
print(dm_dx_sp2)
Ma = np.append(M_in, np.zeros(grid))
TT = np.append(Ttin, np.zeros(grid))
TS = np.append(Ts(Ttin, M_in), np.zeros(grid))

A1 = 0.5612
A2 = 0.1403
A3 = 0.5612
x1 = 0
x2 = 0.2113
x3 = 1
c2 = 3 * (A2 - A1) / (x2 - x1) ** 2
c3 = -2 * (A2 - A1) / (x2 - x1) ** 3
d1 = A2
d2 = 3 * (A3 - A2) / (x3 - x2) ** 2
d3 = -2 * (A3 - A2) / (x3 - x2) ** 3
def dm_f_q(x, Mx, gamma, f, dq_dx):
    if x <= x2:
        Axx = A1 + c2 * (x - x1) ** 2 + c3 * (x - x1) ** 3
        dAxx = 2 * c2 * (x - x1) + 3 * c3 * (x - x1) ** 2
    elif x > x2:
        Axx = d1 + d2 * (x - x2) ** 2 + d3 * (x - x2) ** 3
        dAxx = 2 * d2 * (x - x2) + 3 * d3 * (x - x2) ** 2
    Dxx = np.sqrt(4 * Axx / np.pi)
    XM = (1 + 0.5 * (gamma - 1) * Mx ** 2)
    TTx = Ttin + 398.284*x
    if np.abs(Mx - 1) < 0.01:
        if ((x < x2) and (Mx < 1)) or ((x > x2) and (Mx > 1)):
            return dm_dx_sp2
        else:
            return dm_dx_sp
    else:
        return Mx * (((1 + gamma * Mx ** 2) * XM / (2 * (1 - Mx ** 2))) * (dq_dx / cp / TTx)
                    + (gamma * XM * Mx ** 2 / (2 * (1 - Mx ** 2))) * (f / Dxx)
                    - ((2 + (gamma - 1) * Mx ** 2) / (2 * (1 - Mx ** 2))) * (dAxx / Axx))


sol_sub = solve_ivp(dm_f_q,[0, location], [M_in], t_eval=X[0:ii], args=(gamma, f, dq_dx))
sol_super = solve_ivp(dm_f_q,[location, 1], [1.0], t_eval=X[ii:-1], args=(gamma, f, dq_dx))

fig = plt.figure()
plt.plot(
    np.hstack((sol_sub.t, sol_super.t)),
    np.hstack((sol_sub.y[0, ::], sol_super.y[0, :])), label='with friction and heat addition')
plt.title("mach number vs x for q = 400kj, f = 0.01 and M_in = 0.2")
plt.legend()
plt.grid()
fig.show()

print(location)