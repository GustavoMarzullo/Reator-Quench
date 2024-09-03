import numpy as np
from reator.quench import reator_quench
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def f(X):
    F0 = np.array([0.21825,0.65475,0.05])*6500 #mol/s
    F0N2 = F0[0]
    L1, L2, L3, Pin, T1, Tin, frac1, frac2, frac3 = X
    Y = [frac1, frac2, frac3]
    L2 = L1 + L2
    L3 = L2 + L3
    Lbed = [L1, L2, L3]
    L, T, P, F, rNH3 = reator_quench(Lbed, Pin, Tin, T1, F0, Y,Ac=7)
    FN2, FH2, FNH3 = F[:, 0], F[:, 1], F[:, 2]
    XN2 = (F0N2- FN2[-1])/F0N2
    return -XN2

bounds = [(1, 20), (1, 20), (1, 20), (155,295),(200, 600), (200, 600), (0.01, 0.99), (0.01, 0.99), (0.01, 0.99)]
cons = ({'type': 'eq', 'fun': lambda x:  x[6] + x[7] + x[8] - 1},
        {'type': 'ineq', 'fun': lambda x:  -(x[0] + x[1] + x[2]) + 30})
X0 = [3, 7, 15, 225, 420, 400, 0.5, 0.4, 0.1]
res = minimize(f, X0, bounds=bounds, constraints=cons, options={'disp': True}, tol=1e-6)

print('\n',res.x)
