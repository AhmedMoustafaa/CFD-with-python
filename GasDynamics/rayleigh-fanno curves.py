import numpy as np
import matplotlib.pyplot as plt

gamma = 1.4
cp = 1.0055
M1 = 3  #reference mach number

#reyleigh
M2 = np.linspace(0.000000001,5,100)
ds_r = np.log((M2/M1)**2 * ((gamma*M1**2+1)/(1+gamma*M2**2))**((gamma+1)/gamma))
h2_r = ((1+gamma*M1**2)/(1+gamma*M2**2))**2 * (M2/M1)**2

#fanno
ds_f = np.log((M2 / M1)**((gamma - 1) / gamma) * ((1 + (gamma - 1) / 2 * M1**2) / (1 + (gamma - 1) / 2 * M2**2))**((gamma + 1) / (2 * gamma)))
h2_f = (1+0.5*(gamma-1)*M1**2)/(1+0.5*(gamma-1)*M2**2)



fig, ax = plt.subplots()

ax.plot(ds_r, h2_r, label='Rayleigh Line')
ax.plot(ds_f, h2_f, label='Fanno Line')
ax.set_title(r'Rayleigh And Fanno lines ($\gamma = {gamma}, M_1 = {M1}$)'.format(gamma=gamma,M1=M1))
ax.set_xlabel(r'$\Delta s$')
ax.set_ylabel(r'$H_2$')
ax.grid(True)
ax.set_xlim([-1, 1])
ax.legend()
plt.savefig('fanno-line.png')
plt.tight_layout()
plt.show()
#plt.close()

fig2 = plt.figure()
plt.plot(ds_r,M2)
plt.plot(ds_f,M2)
plt.grid(True)
plt.title(r'Rayleigh And Fanno lines ($\gamma = {gamma}, M_1 = {M1}$)'.format(gamma=gamma,M1=M1))
plt.xlabel(r'$\Delta s$')
plt.ylabel(r'$M$')
plt.xlim([-1, 1])
plt.tight_layout()
plt.savefig('m-s.png')
plt.show()
