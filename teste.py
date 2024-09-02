import numpy as np
from reator.quench import reator_quench
import matplotlib.pyplot as plt

F0 = np.array([0.21825,0.65475,0.05])*6500 #mol/s
Lbed = [3, 10, 50]
Pin = 155 #atm
T1 = 420 #°C
Tin = 390 #ºC
Fin = F0
Y = [0.6, 0.2, 0.2]

####FAZER A PRESSÃO COMO EQUAÇÃO DIFERENCIAL####

L, T, P, F, rNH3 = reator_quench(Lbed, Pin, Tin, T1, Fin, Y, Ac=4)
FN2, FH2, FNH3 = F[:, 0], F[:, 1], F[:, 2]

#https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
fig, axs = plt.subplots(2, 3)
fig.suptitle('Reator químico quench de leito fixo')
axs[0, 0].plot(L, FN2, color='blue')
axs[0, 0].set_title("Vazão de nitrogênio")
axs[0, 0].set(xlabel='Comprimento (m)', ylabel='Vazão (mol/s)')

axs[1, 0].plot(L, FNH3, color='green')
axs[1, 0].set_title("Vazão de amônia")
axs[1, 0].set(xlabel='Comprimento (m)', ylabel='Vazão (mol/s)')

axs[0, 1].plot(L, FH2, color='gray')
axs[0, 1].set_title("Vazão de hidrogênio")
axs[0, 1].set(xlabel='Comprimento (m)', ylabel='Vazão (mol/s)')

axs[1, 1].plot(L, T, color='red')
axs[1, 1].set_title("Temperatura")
axs[1, 1].set(xlabel='Comprimento (m)', ylabel='Temperatura (ºC)')

axs[0, 2].plot(L, P, color='purple')
axs[0, 2].set_title("Pressão")
axs[0, 2].set(xlabel='Comprimento (m)', ylabel='Pressão (atm)')

axs[1, 2].plot(L, rNH3, color='black')
axs[1, 2].set_title("Taxa de reação")
axs[1, 2].set(xlabel='Comprimento (m)', ylabel='rNH3 (mol/(m³.s))')

fig.tight_layout()

plt.show()