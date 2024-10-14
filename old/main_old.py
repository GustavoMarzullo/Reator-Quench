import numpy as np
from Reator.reator_old import reator_quench
import matplotlib.pyplot as plt

F0 = np.array([0.21825,0.65475,0.05])*6500 #mol/s
L1, L2, L3, Pin, T1, Tin, frac1, frac2, frac3 = [2, 5, 20, 180, 420, 400, 0.4, 0.3, 0.2]
Y = [frac1, frac2, frac3]
Lbed = [L1, L2, L3]
L, T, P, F, rNH3 = reator_quench(Lbed, Pin, Tin, T1, F0, Y, Ac=7,method='LSODA')
FN2, FH2, FNH3 = F[:, 0], F[:, 1], F[:, 2]

F0N2 = F0[0]
XN2 = (F0N2- FN2[-1])/F0N2

#https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
fig, axs = plt.subplots(2, 3)
fig.suptitle(f'Reator químico quench de leito fixo. Conversão ={XN2:.2%}')
axs[0, 0].plot(L, FN2, color='black')
axs[0, 0].set_title("Vazão de nitrogênio")
#axs[0, 0].set_ylim(ymin=0)
axs[0, 0].set(xlabel='Comprimento (m)', ylabel='Vazão (mol/s)')

axs[1, 0].plot(L, FNH3, color='black')
axs[1, 0].set_title("Vazão de amônia")
axs[1, 0].set(xlabel='Comprimento (m)', ylabel='Vazão (mol/s)')

axs[0, 1].plot(L, FH2, color='black')
axs[0, 1].set_title("Vazão de hidrogênio")
axs[0, 1].set(xlabel='Comprimento (m)', ylabel='Vazão (mol/s)')

axs[1, 1].plot(L, T, color='black')
axs[1, 1].set_title("Temperatura")
axs[1, 1].set(xlabel='Comprimento (m)', ylabel='Temperatura (ºC)')

axs[0, 2].plot(L, P, color='black')
axs[0, 2].set_title("Pressão")
axs[0, 2].set(xlabel='Comprimento (m)', ylabel='Pressão (atm)')

axs[1, 2].plot(L, rNH3, color='black')
axs[1, 2].set_title("Taxa de reação")
axs[1, 2].set(xlabel='Comprimento (m)', ylabel='rNH3 (mol/(m³.s))')

fig.tight_layout()

plt.show()