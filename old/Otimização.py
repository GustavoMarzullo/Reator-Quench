import numpy as np
from reator.quench import reator_quench
import matplotlib.pyplot as plt
from random import uniform
from scipy.optimize import minimize
import warnings
from time import time
warnings.filterwarnings("ignore")

def f(X):
    F0 = np.array([0.21825,0.65475,0.05])*6500 #mol/s
    F0N2 = F0[0]
    L1, L2, L3,L4, Pin, T1, Tin, frac1, frac2, frac3, frac4 = X
    Y = [frac1, frac2, frac3, frac4]
    L2 = L1 + L2
    L3 = L2 + L3
    L4 = L3 + L4
    Lbed = [L1, L2, L3, L4]
    L, T, P, F, rNH3 = reator_quench(Lbed, Pin, Tin, T1, F0, Y,Ac=7,method='LSODA')
    FN2, FH2, FNH3 = F[:, 0], F[:, 1], F[:, 2]
    XN2 = (F0N2- FN2[-1])/F0N2
    return -XN2

bounds = [(1, 20), (1, 20), (1, 20), (1, 20) ,(155,295),(400, 500), (350, 400), (0.01, 0.99), (0.01, 0.99), (0.01, 0.99), (0.01, 0.99)]
cons = ({'type': 'eq', 'fun': lambda x:  x[7] + x[8] + x[9] + x[10] - 1},
        {'type': 'ineq', 'fun': lambda x:  -(x[0] + x[1] + x[2]) + 40})
#X0 = [1, 1, 20, 225, 420, 400, 0.5, 0.4, 0.1]

def optimize():
    X = []
    for i in range(10):
        print(f"Tentativa {i+1}")
        X0 = [uniform(i[0],i[1]) for i in bounds]

        X0[1] = uniform(0, 40-X0[0])
        X0[2] = uniform(0, 40-X0[1]-X0[0])
        X0[3]= 40 - X0[2] - X0[1] -X0[0]

        X0[-3] = uniform(0, 1-X0[-4])
        X0[-2] = uniform(0, 1-X0[-3]-X0[-4])
        X0[-1]= 1 - X0[-2] - X0[-3] -X0[-4]
        print("Chute inicial:")
        print(X0)
        try:
            res = minimize(f, X0, bounds=bounds, constraints=cons, options={'disp': True})
            print(res.x)
            X.append([i+1,f(res.x), res.x])
            print("\n")
        except:
            X.append([i+1,0, 0])
            print("Erro\n")
            continue
    print("Melhor resultado:")
    print(min(X, key=lambda x: x[1]))
    print("\n")
optimize()
