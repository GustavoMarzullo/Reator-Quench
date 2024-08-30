import numpy as np
from reator.quench import reator_quench
import matplotlib.pyplot as plt

F0 = np.array([0.21825,0.65475,0.05])*6500 #mol/s
P0 = 155 #atm
T0 = 430 #ºC
X = 0
L = 0
Lbed = [10,20,30]
Pin = 155 #atm
T1 = 400 #°C
Tin = 250 #ºC
Fin = F0
Y = [0.8,0.15,0.05]


L, T, F = reator_quench(Lbed, Pin, Tin, T1, Fin, Y)
FN2, FH2, FNH3 = F[:, 0], F[:, 1], F[:, 2]

#https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(L, FN2, color='blue')
axs[0, 0].set_title("Vazão de nitrogênio")
axs[1, 0].plot(L, FNH3, color='green')
axs[1, 0].set_title("Vazão de amônia")
axs[0, 1].plot(L, FH2, color='gray')
axs[0, 1].set_title("Vazão de hidrogênio")
axs[1, 1].plot(L, T, color='red')
axs[1, 1].set_title("Temperatura")
fig.tight_layout()

plt.show()