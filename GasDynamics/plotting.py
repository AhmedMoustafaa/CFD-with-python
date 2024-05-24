import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

X = np.linspace(0,1,10000+1)

df1 = pd.read_csv('arrays1.csv')
df2 = pd.read_csv('arrays2.csv')
df3 = pd.read_csv('arrays3.csv')

M1 = df1['Mach1'].to_numpy()
M2 = df2['Mach2'].to_numpy()
M3 = df3['Mach3'].to_numpy()

V1 = df1['V1'].to_numpy()
V2 = df2['V2'].to_numpy()
V3 = df3['V3'].to_numpy()

TT1 = df1['TT1'].to_numpy()
TT2 = df2['TT2'].to_numpy()
TT3 = df3['TT3'].to_numpy()

PT1 = df1['PT1'].to_numpy()
PT2 = df2['PT2'].to_numpy()
PT3 = df3['PT3'].to_numpy()

TS1 = df1['Ts1'].to_numpy()
TS2 = df2['Ts2'].to_numpy()
TS3 = df3['Ts3'].to_numpy()

PS1 = df1['PS1'].to_numpy()
PS2 = df2['PS2'].to_numpy()
PS3 = df3['PS3'].to_numpy()

fig1 = plt.figure()
plt.plot(X,M1,label='area change + friction + Heat addition')
plt.plot(X,M2,label='area change + friction')
plt.plot(X,M3,label='area change')
plt.grid(True)
plt.title(r'Mach number ($\gamma = {gamma}, M_1 = {M1}$)')
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

fig3 = plt.figure()
plt.plot(X,PT1,label='area change + friction + Heat addition')
plt.plot(X,PT2,label='area change + friction')
plt.plot(X,PT3,label='area change')
plt.grid(True)
plt.title(r'Pressure')
plt.legend()
fig3.savefig('pressure')

fig4 = plt.figure()
plt.plot(X,TS1,label='area change + friction + Heat addition')
plt.plot(X,TS2,label='area change + friction')
plt.plot(X,TS3,label='area change')
plt.grid(True)
plt.title('static Temperature')
plt.legend()
fig4.savefig('static temperature')

fig5 = plt.figure()
plt.plot(X,PS1,label='area change + friction + Heat addition')
plt.plot(X,PS2,label='area change + friction')
plt.plot(X,PS3,label='area change')
plt.grid(True)
plt.title('static pressure')
plt.legend()
fig5.savefig('static_pressure')

fig6 = plt.figure()
plt.plot(X,V1,label='area change + friction + Heat addition')
plt.plot(X,V2,label='area change + friction')
plt.plot(X,V3,label='area change')
plt.grid(True)
plt.title(r'Velocity')
plt.legend()
fig6.savefig('velocity')