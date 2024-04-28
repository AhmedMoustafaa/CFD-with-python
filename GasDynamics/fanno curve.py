import numpy as np
from matplotlib import pyplot as plt

gamma = 1.4
M1 = 5.0
Tt_T1 = 1 + 0.5*(gamma-1)*M1**2

max_temp_ratio = (1+gamma*M1**2)**2 / (4*gamma*M1**2)

T_T1 = np.linspace(1.0, max_temp_ratio, 300)
delta_s_cp_upper = (
    np.log(T_T1) - (gamma-1)*np.log(
        (1+gamma*M1**2 + np.sqrt((1+gamma*M1**2)**2 - 4*gamma*T_T1*M1**2))/2
        ) / gamma
    )
delta_s_cp_lower = (
    np.log(T_T1) - (gamma-1)*np.log(
        (1+gamma*M1**2 - np.sqrt((1+gamma*M1**2)**2 - 4*gamma*T_T1*M1**2))/2
        ) / gamma
    )

delta_s_cp = (1/gamma)*np.log(T_T1) + (gamma-1)*np.log((Tt_T1 - T_T1)/(Tt_T1 - 1)) / (2*gamma)


fig = plt.subplots()
plt.plot(np.hstack((delta_s_cp_lower, delta_s_cp_upper[::-1])),
    np.hstack((T_T1, T_T1[::-1])))
plt.plot(delta_s_cp,T_T1)
plt.title(r'Fanno line ($\gamma = {gamma}$)'.format(gamma=gamma))
plt.xlabel(r'$\Delta s / c_p$')
plt.ylabel(r'$T / T_1$')
plt.grid(True)
#plt.savefig('fanno-line.png')
plt.tight_layout()
plt.show()
#plt.close()