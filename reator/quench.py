import numpy as np
import sys
sys.path.append('reator')
from reator.conversao import dX_dL
from reator.temperatura import dT_dL
from pressao import ergun
from scipy.integrate import solve_ivp
import numpy as np

def calc_ODE(
        L: float, 
        Y: np.ndarray, 
        F0: np.ndarray, 
        P0: float) -> np.ndarray:
    """
    Função para calcular a conversão e a temperatura em relação ao comprimento do reator.
    
    Parameters
    ----------
    L 
        Comprimento do reator [m].
    Y 
        Vetor com a conversão e temperatura iniciais [adimensional, ºC].
    F0 
        Vazão de reagentes na entrada do reator (N2, H2, NH3) [mol/s].
    P0 
        Pressão na entrada do reator [atm].
    
    Retorna
    -------
    np.array
        Conversão e temperatura em relação ao comprimento do reator. [adimensional, ºC].
    """
    X, T = Y
    derivada_conversao = dX_dL(X, L, T, F0, P0)[0]
    derivada_temperatura = dT_dL(X, L, T, F0, P0)
    return np.array([derivada_conversao, derivada_temperatura])

def reator_quench(Lbed, Pin, Tin, T1, Fin, Y, number_of_points=31):
    L, T, F, P = np.array([]), np.array([]), np.empty([3]), np.array([])
    assert len(Lbed) == len(Y)
    for i in range(len(Lbed)):
        if i==0:
            #resolvendo a EDO
            T0 = T1
            Y0 = [0, T0]
            L_eval = np.linspace(0, Lbed[i], number_of_points)
            F0 = Fin*Y[i]
            P0 = Pin
            sol = solve_ivp(calc_ODE, [0, Lbed[i]], Y0, t_eval=L_eval, args=(F0, P0))
            Xsol, Tsol = sol.y[0], sol.y[1]
            #calculando as composições
            F0N2, F0H2, F0NH3 = F0
            FN2 = F0N2*(1-Xsol)
            FH2 = F0H2 - 3*(F0N2-FN2)
            FNH3 = F0NH3 + 2*(F0N2-FN2)
            Fsol = np.array(list(zip(FN2, FH2, FNH3)))
            #criando os vetores 
            L, T = np.append(L, L_eval), np.append(T, Tsol)
            for _ in range(len(Fsol)):
                F = np.vstack([F, Fsol[_]])
            #calculando a pressão
            for j in range(len(L_eval)):
                Pcalc = ergun(L_eval[j],P0,F[j])
                P = np.append(P, Pcalc)
        else:
            T0 = (T[-1]*sum(F[-1]) + Tin*sum(Fin)*Y[i])/(sum(F[-1]) + sum(Fin)*Y[i]) #isso pode ser melhorado!!
            Y0 = [0, T0]
            L_eval = np.linspace(L[-1], Lbed[i], number_of_points)
            F0 = Fin*Y[i] + F[-1] 
            P0 = P[-1]
            sol = solve_ivp(calc_ODE, [L[-1], Lbed[i]], Y0, t_eval=L_eval, args=(F0, P0))
            Xsol, Tsol = sol.y[0], sol.y[1]
            #calculando as composições
            F0N2, F0H2, F0NH3 = F0
            FN2 = F0N2*(1-Xsol)
            FH2 = F0H2 - 3*(F0N2-FN2)
            FNH3 = F0NH3 + 2*(F0N2-FN2)
            Fsol = np.array(list(zip(FN2, FH2, FNH3)))
            #criando os vetores 
            L, T = np.append(L, L_eval), np.append(T, Tsol)
            for _ in range(len(Fsol)):
                F = np.vstack([F, Fsol[_]])
            #calculando a pressão
            for j in range(len(L_eval)):
                Pcalc = ergun(L_eval[j],P0,F[j])
                P = np.append(P, Pcalc)

    return L, T, F[1:len(F)]

            
