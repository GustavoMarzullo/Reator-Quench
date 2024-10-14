import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('Reator')
from Reator.conversao import dX_dL
from Reator.temperatura import dT_dL
from pressao import dP_dL
from scipy.integrate import solve_ivp
from dataclasses import dataclass

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

@dataclass
class Reagente:
    FN2: float
    FH2: float
    FNH3: float
    P: float
    T: float

    def F_total(self):
        '''Returns the sum of the reagent flows'''
        return self.FN2 + self.FH2 + self.FNH3
    
    def F(self):
        return np.array([self.FN2, self.FH2, self.FNH3])

@dataclass
class Leito:
    """ 
    L 
        Comprimento do leito [m].
    T
        Temperatura de entrada no leito [°C].
    φ 
        Porosidade do leito fixo [adimensional].
    Ac 
        Área da seção transversal do leito [m²].
    Dp
        Diâmetro das partilhas do leito fixo [mm]."""
    
    L: float
    T:float | None = None
    φ:float=0.4
    Ac:float=7
    Dp:float=2

@dataclass
class Reator:
    """
    Leitos: [Leito]
        Leitos do reator

    Reagente: Reagente
        Carga do reator

    Y: np.ndarray
        Fração de carga direcionada a cada leito
    """
    Leitos: [Leito]  # type: ignore
    Reagente: Reagente
    Y: np.ndarray

    def solve(self,  number_of_points=11, metodo='LSODA'):
        """
        Resolve o sistema de equações diferenciais para o reator químico.

        Parameters
        ----------
        number_of_points : int, optional
            Número de pontos para a resolução da EDO. O padrão é 11.
        method : str, optional
            Método de resolução da EDO. O padrão é 'LSODA'.
            
        Retorna
        -------
        L : numpy.ndarray
            Comprimento do reator [m]
        T : numpy.ndarray
            Temperatura do reator [°C]
        P : numpy.ndarray
            Pressão do reator [atm]
        F : numpy.ndarray
            Composição do reator (N2, H2, NH3) [mol/s]
        rNH3 : numpy.ndarray
            Taxa de reação [mol/s]"""
    
        assert len(self.Leitos) == len(self.Y), "O número de leitos deve ser o mesmo que o de frações de carga."
        assert np.isclose(sum(self.Y), 1) , "A soma das frações de carga deve ser igual a 1."

        L, T, F, P, rNH3, XN2 = np.array([]), np.array([]), np.empty([3]), np.array([]), np.array([]), np.array([])

        Lin = 0.0
        for i in range(len(self.Y)):
            #calculando os valores iniciais do sistema de EDO
            Tin = self.Reagente.T if i ==0 else self.Leitos[i].T
            if i!=0:
                Fin = self.Reagente.F_total()*self.Y[i]
                T0 = (T[-1]*sum(F[-1]) + Tin*Fin)/(sum(F[-1]) + Fin)
            else:
                T0 = Tin
            P0 = self.Reagente.P if i==0 else P[-1]
            F0 = self.Reagente.F()*self.Y[i] if i==0 else self.Reagente.F()*self.Y[i] + F[-1]
            L_eval = np.linspace(0, self.Leitos[i].L, number_of_points)
            Y0 = [0, Tin, P0]
            φ, Ac, Dp = self.Leitos[i].φ, self.Leitos[i].Ac, self.Leitos[i].Dp

            #resolvendo o sistema de EDO
            sol = solve_ivp(calc_ODE, [0, self.Leitos[i].L], Y0, t_eval=L_eval, args=(F0, φ, Ac, Dp), method=metodo)
            Xsol, Tsol, Psol = sol.y[0], sol.y[1], sol.y[2]

            #calculando as composições
            F0N2, F0H2, F0NH3 = F0
            FN2 = F0N2*(1-Xsol)
            FH2 = F0H2 - 3*(F0N2-FN2)
            FNH3 = F0NH3 + 2*(F0N2-FN2)

            #salvando os resultados
            Fsol = np.array(list(zip(FN2, FH2, FNH3)))
            for _ in range(len(Fsol)):
                F = np.vstack([F, Fsol[_]])
            for j in range(len(L_eval)):
                derivada, _rNH3, _F, Pcalc, η = dX_dL(Xsol[j], Psol[j], Tsol[j], F0, Ac)
                rNH3 = np.append(rNH3, _rNH3)
            L, T, P = np.append(L, (L_eval+Lin)), np.append(T, Tsol), np.append(P, Psol)

            Lin = Lin + self.Leitos[i].L

        resultado = Resultados(L, T, P, F[1:len(F)], rNH3, self.Reagente)

        return resultado


@dataclass
class Resultados:
    L: np.ndarray
    T: np.ndarray
    P: np.ndarray
    F: np.ndarray
    rNH3: np.ndarray
    Reagente: Reagente
    
    def FN2(self):
        return self.F[:, 0]

    def FH2(self):
        return self.F[:, 1]
    
    def FNH3(self):
        return self.F[:, 2]
    
    def XN2(self):
        F0N2 = self.Reagente.FN2
        return (F0N2 - self.FN2()[-1]) / F0N2
    
    def plot(self):
        L = self.L
        FN2, FH2, FNH3, T, P, XN2, rNH3 = self.FN2(), self.FH2(), self.FNH3(), self.T, self.P, self.XN2(), self.rNH3
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
    



