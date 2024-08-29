import numpy as np
import sys
sys.path.append('reator')
from conversao import dX_dL

def cP_mix(
        P:float, 
        T: float) -> float:
    """
    Calcula o cP de uma mistura de gases para a síntese de amônia.
    Fonte: Suhan et. atl (2022)

    Parameters
    ----------
    P 
        Pressão [atm].
    T 
        Temperatura [K].

    Retorna
    -------
    cP : float
        cP de uma mistura de gases para a síntese de amônia [J/(mol.ºC)]
    """
    P = P * 101.325 #atm -> kPa

    cP = 35.31 + 0.02*T  + 0.00000694*T**2 - 0.0056*P + 0.000014*P*T
    
    return cP

def calor_de_reacao(P, T): 
    """
    Calcula o calor de reação da reação de síntese de amônia.
    Fonte: Gunorubon (2014)

    Parameters
    ----------
    P 
        Pressão [atm].
    T 
        Temperatura [K].

    Returns
    -------
    ΔHr : float
        Calor de reação da reação de síntese de amônia [J/mol]"""
    
    ΔHr = 4.184*(-(0.54526 + 840.609/T + 459.374e6/T**3)*P - 5.34685*T - 0.2525e-3*T**2 + 1069197e-6*T**3 - 9157.09)

    return -41.14e3 #ΔHr 
      

def dT_dL(
        X:float,
        L: float,
        T:float,
        F0:np.ndarray,
        P0:float,
        φ:float=0.4,
        Ac:float=7)-> float:
        """
        Função para calcular a derivada da temperatura no reator.
        Fonte: Jorqueira et. al (2018) 

        Parameters
        ----------
        X
            Conversão de nitrogênio [admensional].
        L
            Comprimento do reator [m].
        T
            Temperatura [ºC].
        F0
            Vazão de reagentes na entrada do reator (N2, H2, NH3) [mol/s].
        P0
            Pressão na entrada do reator [atm].
        φ
            Porosidade do leito fixo [adimensional].
        Ac
            Área da seção transversal do reator [m²].
        Returns
        -------
        Derivada da temperatura no reator [ºC/m]."""

        _dX_dL, rNH3, F, P, η = dX_dL(X, L, T, F0, P0, φ, Ac)

        T = T + 273.15 #ºC -> K
        cP = cP_mix(P, T)
        ΔHr = calor_de_reacao(P, T)

        dT_dL = (η*(-ΔHr)*Ac*rNH3)/(sum(F)*cP)

        return dT_dL

