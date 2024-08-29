import numpy as np
#Cp = 3200 J/(kg.K)

def coeficiente_atividade(
        P: float,
        T: float) -> np.ndarray:
    """
    Função para calcular o coeficiente de atividade do nitrogênio, hidrogênio e amônia.
    Fonte: Dyson (1968)

    Parameters
    ----------
    P
        Pressão [atm].
    T
        Temperatura [K].
    
    Returns
    -------
    Coeficiente de atividade do nitrogênio, hidrogênio e amônia, respectivamente [adimensional]
    """
    γN2 = 0.93431737 + 0.310180e-3*T + 0.295896e-3*P - 0.270727e-6*T**2 + 0.4775207e-6*P**2
    γH2 = np.exp(np.exp(-3.8402*T**0.125 +0.541)*P - np.exp(-0.1263*T**0.5 -15.980)*P**2 + 300*(np.exp(-0.011901*T-5.491)*(np.exp(-P/300)-1)))
    γNH3 = 0.1438996 +  0.2028538e-2*T - 0.4487672e-3*P - 0.1142945e-5*T**2 + 0.2761216e-6*P**2


    return np.array([γN2, γH2, γNH3])

def constante_velocidade(T: float) -> float:
    """
    Função para calcular a constante de velocidade da formação da amônia.

        Parameters
    ----------
    T
        Temperatura [K].
    
    Returns
    -------
    Constante de velocidade da formação da amônia [adimensional]
    """
    k = 2.457e11*np.exp(-40765/(1.987*T))
    return k

def constante_equilibrio_amonia(T:float) -> float:
    """
    Função para calcular a constante de equilibrio da amônia.
    Fonte: Dyson (1968)

        Parameters
    ----------
    T
        Temperatura [K].
    
    Returns
    -------
    Constante de equilíbrio da amônia [adimensional]
    """
    Ka = 10**(-2.691122*np.log10(T) - 5.519265e-5*T + 1.848863e-7*T**2 + 2001.6*T**(-1) + 2.6899)
    return Ka

###constantes do fator de efetividade do nitrogênio
B0 = [-17.539096, -8.2125534, -4.6757259]
B1 = [0.07697849, 0.03774149, 0.02354872]
B2 = [6.900548, 6.190112, 4.687353]
B3 = [-1.082790e-4, -5.354571e-5, -3.463308e-5]
B4 = [-26.424699, -20.86963, -11.28031]
B5 = [4.927648e-8, 2.379142e-8, 1.540881e-8]
B6 = [38.93727,	27.88403, 10.46627]
Pressoes_B = [150, 225, 300]
###

def constantes_fator_efetividade_N2(P:float) -> np.ndarray:
    """
    Função para calcular as constantes do fator de efetividade do nitrogênio.
    Fonte: Dyson (1968)

        Parameters
    ----------
    P
        Pressão [atm].
    
    Returns
    -------
    Constantes do fator de efetividade do nitrogênio [adimensional]
    """
    b0 = np.interp(P, Pressoes_B, B0)
    b1 = np.interp(P, Pressoes_B, B1)
    b2 = np.interp(P, Pressoes_B, B2)
    b3 = np.interp(P, Pressoes_B, B3)
    b4 = np.interp(P, Pressoes_B, B4)
    b5 = np.interp(P, Pressoes_B, B5)
    b6 = np.interp(P, Pressoes_B, B6)

    B = np.array([b0, b1, b2, b3, b4, b5, b6])

    return B

def fator_efetividade_N2(
        P: float,
        T: float,
        X: float) -> float:
    """
    Função para calcular o fator de efetividade do nitrogênio.
    Fonte: Dyson (1968)

    Parameters
    ----------
    P
        Pressão [atm].
    T
        Temperatura [K].
    X
        Conversão de amônia [adimensional].
    
    Returns
    -------
    Fator de efetividade do nitrogênio. [adimensional]
    """
    B = constantes_fator_efetividade_N2(P)
    η =  B[0] + B[1]*T + B[2]*X + B[3]*T**2 + B[4]*X**2 + B[5]*T**3 + B[6]*X**3

    return η

def dP_dL( 
    P:float,
    T: float,
    X:float,
    P0:float,
    T0:float,
    F0:np.ndarray,
    Ac: float = 7,
    φ:float=0.4,
    μ:float= 0.0275,
    ρ:float = 23,
    ρc:float = 1816.5,
    Dp:float=7,
    δ = -0.5) -> float:
    """
    Função para calcular a pressão em um reator químico de leito fixo.

    Parameters
    ----------
    P
        Pressão  do reator [atm].
    T
        Temperatura do reator [ºC].
    X
        Conversão de N2 para amônia [adimensional].
    P0
        Pressão na entrada do reator [atm].
    T0
        Temperatura na entrada do reator [ºC]
    F0
        Vazão de reagentes na entrada do reator (N2, H2, NH3) [mol/(s)].
    Ac
        Área da seção transversal do reator [m²].
    φ
        Porosidade do leito fixo [adimensional].
    μ
        Viscosidade do gás [cP].
    ρ
        Densidade do gás [kg/m³].
    ρc
        Densidade do catalisador [kg/m³].     
    Dp
        Diâmetro das partículas do leito fixo [mm].
    δ
        Variação no número total de mols por mol de N2  reagido

    Returns
    -------
    Derivada da pressão na saída do reator [atm/m].
    """
    T = T+273.15 #°C -> K
    T0 = T0+273.15 #°C -> K
    μ = μ*0.001 #cP -> Pa.s

    yN20 = F0[0]/sum(F0)
    ε = δ*yN20
    Dp = Dp*1e-3 #mm -> m
    G = sum(F0) * 4.25e-3 #velocidade mássica superficial (vazão na entrada * massa molar média na entrada)


    β0 =  G*(1-φ)/(ρ*Dp*φ**3) * ((150*(1-φ)*μ)/Dp + 1.75*G)
    α = (2*β0)/(Ac*ρc*(1-φ)*P0)

    dP_dL = ((1-φ)*Ac*ρc) * -α/2 * T/T0 * P0/(P-P0) * (1+ ε*X)

    return dP_dL

def dX_dL(
        P:float,
        T:float,
        X:float,
        F0:np.ndarray,
        φ:float=0.4,
        Ac:float=7,         
        ρc:float=1816.5):    

    """
    Função para calcular a derivada da conversão de N2 para amônia.

    Parameters
    ----------
    P
        Pressão [atm].
    T
        Temperatura [K].
    X
        Conversão de nitrogênio [admensional].
    F0
        Vazão de reagentes na entrada do reator (N2, H2, NH3) [mol/s].
    φ
        Porosidade do leito fixo [adimensional].
    Ac
        Área da seção transversal do reator [m²].
    ρc
        Densidade do catalisador [kg/m³].
    Returns
    -------
    Derivada da conversão de nitrogênio em relação ao comprimento do reator [1/s].
    """
    T = T+273.15 #°C -> K
    Ka = constante_equilibrio_amonia(T)
    k = constante_velocidade(T)
    aN2, aH2, aNH3 = coeficiente_atividade(P, T) * F0
    η = fator_efetividade_N2(P, T, X)
    rNH3 = 2*k*(Ka**2 * aN2*(aH2**3/aNH3**2)**0.5 - (aNH3**2/aH2**3)**0.5) * η
    rN2 = -1/2*rNH3

    dX_dL = ((1-φ)*Ac*ρc)/sum(F0) * rN2
    return dX_dL





