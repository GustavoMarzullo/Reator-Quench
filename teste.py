import numpy as np
from reator.pressao import ergun
from reator.conversao import dX_dL, constante_equilibrio_amonia
from reator.temperatura import dT_dL

F0 = np.array([0.21825,0.65475,0.05])*6500 #mol/s
P0 = 155 #atm
T0 = 430 #ºC
X = 0
L = 0

#print(dX_dL(X, L, T0, F0, P0))
#print(dT_dL(X, L, T0, F0, P0))
print(constante_equilibrio_amonia(500))