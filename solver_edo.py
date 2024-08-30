from reator.conversao import dX_dL
from reator.temperatura import dT_dL
from reator.pressao import ergun
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def calc_ODE(L, Y, F0, P0):
    X, T = Y
    derivada_conversao = dX_dL(X, L, T, F0, P0)[0]
    derivada_temperatura = dT_dL(X, L, T, F0, P0)
    return np.array([derivada_conversao, derivada_temperatura])

L0, Lf = 0, 30 #m
F0 = np.array([0.2173,0.6608,0.0380])*23318e3/3600 #mol/s
X0 = 0
P0 = 155 #atm
T0 = 400 #ºC

L_eval = np.linspace(L0, Lf,101)
Y0 = np.array([X0, T0]) #conversão, temperatura

sol = solve_ivp(calc_ODE, [L0, Lf], Y0, t_eval=L_eval, args=(F0, P0))

Xout = sol.y[0]
Tout = sol.y[1]
print(Xout[-1], Tout[-1])
fig, ax = plt.subplots()
ax.plot(sol.t, Xout, color='black', linewidth=1)
ax1 = ax.twinx()
ax1.plot(sol.t, Tout, color='red', linewidth=1)
ax.set_xlabel("Comprimento (m)")
ax.set_ylabel("Conversão (admensional)", color='black')
ax1.set_ylabel("Temperatura (ºC)", color='red')
plt.show()


