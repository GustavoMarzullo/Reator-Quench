import numpy as np
from reator.quench import reator_quench
import matplotlib.pyplot as plt

F0 = np.array([0.21825,0.65475,0.05])*6500 #mol/s
Lbed = [3, 6, 30]
Pin = 255 #atm
T1 = 420 #°C
Tin = 400 #ºC
Fin = F0
Y = [1/2, 1/4, 1/4]

####FAZER A PRESSÃO COMO EQUAÇÃO DIFERENCIAL####

L, T, P, F, rNH3 = reator_quench(Lbed, Pin, Tin, T1, Fin, Y, Ac=4)
FN2, FH2, FNH3 = F[:, 0], F[:, 1], F[:, 2]

#https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
fig, axs = plt.subplots(2, 3)
axs[0, 0].plot(L, FN2, color='blue')
axs[0, 0].set_title("Vazão de nitrogênio")
axs[1, 0].plot(L, FNH3, color='green')
axs[1, 0].set_title("Vazão de amônia")
axs[0, 1].plot(L, FH2, color='gray')
axs[0, 1].set_title("Vazão de hidrogênio")
axs[1, 1].plot(L, T, color='red')
axs[1, 1].set_title("Temperatura")
axs[0, 2].plot(L, P, color='purple')
axs[0, 2].set_title("Pressão")
axs[1, 2].plot(L, rNH3, color='black')
axs[1, 2].set_title("Taxa de reação")
fig.tight_layout()

plt.show()