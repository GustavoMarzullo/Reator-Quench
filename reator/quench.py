import numpy as np
import sys
sys.path.append('reator')
from reator.conversao import dX_dL
from reator.temperatura import dT_dL
from pressao import dP_dL
from scipy.integrate import solve_ivp

def calc_ODE(
        L: float, 
        Y: np.ndarray, 
        F0: np.ndarray, 
        φ:float=0.4,
        Ac:float=7,
        Dp:float=2) -> np.ndarray: 
    """
    Função para calcular a conversão e a temperatura em relação ao comprimento do reator.
    
    Parameters
    ----------
    L 
        Comprimento do reator [m].
    Y 
        Vetor com a conversão, temperatura e pressão iniciais [adimensional, ºC].
    F0 
        Vazão de reagentes na entrada do reator (N2, H2, NH3) [mol/s].
    φ 
        Porosidade do leito fixo [adimensional].
    Ac 
        Área da seção transversal do reator [m²].
    Dp
        Diâmetro das partilhas do leito fixo [mm].
    
    Retorna
    -------
    np.array
        Conversão e temperatura em relação ao comprimento do reator. [adimensional, ºC].
    """
    X, T, P = Y
    derivada_conversao, _rNH3, _F, _P, _η = dX_dL(X, P, T, F0, Ac)
    derivada_temperatura = dT_dL(X, P, T, F0, Ac)
    derivada_pressao = dP_dL(_F, Ac, φ, Dp)

    return np.array([derivada_conversao, derivada_temperatura, derivada_pressao])

def reator_quench(
        Lbed:np.ndarray, 
        Pin:float, 
        Tin:float, 
        T1:float, 
        Fin:np.ndarray, 
        Y:np.ndarray,
        φ:float=0.4,
        Ac:float=7,
        Dp = 2,
        number_of_points=31) -> tuple:
    """
    Resolve a equação diferencial para um reator químico de leito fixo.

    Parameters
    ----------
    Lbed : 
        Comprimento de cada leito do reator [m].
    Pin : 
        Pressão na entrada do reator [atm].
    Tin : 
        Temperatura na entrada do reator [ºC].
    T1 : 
        Temperatura de entrada no primeiro leito [ºC].
    Fin : 
        Vazão de reagentes na entrada do reator (N2, H2, NH3) [mol/s].
    Y : 
        Fração da carga que é direcionada a cada leito [adimensional].
    φ 
        Porosidade do leito fixo [adimensional].
    Ac 
        Área da seção transversal do reator [m²].
    Dp
        Diâmetro das partilhas do leito fixo [mm].
    number_of_points :
        Número de pontos para a resolução da EDO. O padrão é 31.

    Retorna
    -------
    L, T, P, F :
        Comprimento do reator [m], temperatura [ºC], pressão [atm], composição do reator (N2, H2, NH3) [mol/s] e taxa de reação
    """
    L, T, F, P, rNH3, XN2 = np.array([]), np.array([]), np.empty([3]), np.array([]), np.array([]), np.array([])
    assert len(Lbed) == len(Y)
    for i in range(len(Lbed)):
        if i==0:
            #resolvendo a EDO
            T0 = T1
            P0 = Pin
            Y0 = [0, T0, P0]
            L_eval = np.linspace(0, Lbed[i], number_of_points)
            F0 = Fin*Y[i]
            sol = solve_ivp(calc_ODE, [0, Lbed[i]], Y0, t_eval=L_eval, args=(F0, φ, Ac, Dp))
            Xsol, Tsol, Psol = sol.y[0], sol.y[1], sol.y[2]
            #calculando as composições
            F0N2, F0H2, F0NH3 = F0
            FN2 = F0N2*(1-Xsol)
            FH2 = F0H2 - 3*(F0N2-FN2)
            FNH3 = F0NH3 + 2*(F0N2-FN2)
            Fsol = np.array(list(zip(FN2, FH2, FNH3)))
            for _ in range(len(Fsol)):
                F = np.vstack([F, Fsol[_]])
            for j in range(len(L_eval)): 
                derivada, _rNH3, _F, Pcalc, η = dX_dL(Xsol[j], Psol[j], Tsol[j], F0, Ac)
                rNH3 = np.append(rNH3, _rNH3)
            #criando os vetores 
            L, T, P = np.append(L, L_eval), np.append(T, Tsol), np.append(P, Psol)
        else:
            T0 = (T[-1]*sum(F[-1]) + Tin*sum(Fin)*Y[i])/(sum(F[-1]) + sum(Fin)*Y[i]) #isso pode ser melhorado!!
            L_eval = np.linspace(L[-1], Lbed[i], number_of_points)
            F0 = Fin*Y[i] + F[-1] 
            P0 = P[-1]
            Y0 = [0, T0, P0]
            sol = solve_ivp(calc_ODE, [L[-1], Lbed[i]], Y0, t_eval=L_eval, args=(F0, φ, Ac, Dp))
            Xsol, Tsol, Psol = sol.y[0], sol.y[1], sol.y[2]
            #calculando as composições
            F0N2, F0H2, F0NH3 = F0
            FN2 = F0N2*(1-Xsol)
            FH2 = F0H2 - 3*(F0N2-FN2)
            FNH3 = F0NH3 + 2*(F0N2-FN2)
            Fsol = np.array(list(zip(FN2, FH2, FNH3)))
            for _ in range(len(Fsol)):
                F = np.vstack([F, Fsol[_]])
            for j in range(len(L_eval)):
                derivada, _rNH3, _F, Pcalc, η = dX_dL(Xsol[j], Psol[j], Tsol[j], F0, Ac)
                rNH3 = np.append(rNH3, _rNH3)
            L, T, P = np.append(L, L_eval), np.append(T, Tsol), np.append(P, Psol)

    return L, T, P, F[1:len(F)], rNH3

            
