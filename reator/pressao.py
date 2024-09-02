import numpy as np

def dP_dL( 
    F:np.ndarray,
    Ac: float = 7,
    φ:float=0.4,
    Dp:float=2) -> float:
    """
    Calcula a derivada da pressão no reator químico de leito fixo.
    Fonte: Jorqueira (2018)

    Parameters
    ----------
    F 
        Vazão de reagentes do reator (N2, H2, NH3) [mol/s].
    Ac 
        Área da seção transversal do reator [m²].
    φ 
        Porosidade do leito fixo [adimensional].
    Dp 
        Diâmetro das partículas do leito fixo [mm].

    Returns
    -------
    float
        Derivada da pressão no reator [atm/m].
    """
    μ = 2.254e-5 #Pa.s
    ρ = 23 #kg/m³
    Dp = Dp*1e-3 #mm -> m
    #P0 = P0*101325 #atm -> Pa
    #T = T+273.15 #°C -> K
    #T0 = T0+273.15 #°C -> K

    VazaoVolumetrica =(sum(F)*9.437e-3)/ρ #mol/s -> m³/s
    v = VazaoVolumetrica/Ac #m³/s -> m/s
    G = v*φ #velocidade superficial do gás

    #β0 = (G*(1-φ))/(ρ*Dp*φ**3)*((150*(1-φ)*μ)/Dp +1.75*G)
    #derivada = -β0 #* P0/P * T/T0 * sum(F)/sum(F0)

    derivada = (150*(1-φ)**2*μ*G)/(φ**3*Dp**2) - (1.75*(1-φ)*ρ*G**2)/(φ**3*Dp) #Pa/m
    derivada = derivada/101300 #Pa/m -> atm/m

    return derivada

def ergun( 
    L: float,
    P0:float,
    F:np.ndarray,
    Ac: float = 7,
    φ:float=0.4,
    Dp:float=2) -> float:
    """
    Calcula a pressão em um reator químico de leito fixo.
    Fonte: Suhan et. al (2022)

    Parameters
    ----------
    L 
        Comprimento do reator [m].
    P0 
        Pressão na entrada do reator [atm]. 
    F 
        Vazão de reagentes do reator (N2, H2, NH3) [mol/s].
    Ac 
        Área da seção transversal do reator [m²].
    φ 
        Porosidade do leito fixo [adimensional].
    Dp 
        Diâmetro das partículas do leito fixo [mm].

    Returns
    -------
    float
        Pressão na saída do reator [atm].
    """
    μ = 2.254e-5 #Pa.s
    ρ = 23 #kg/m³
    Dp = Dp*1e-3 #mm -> m
    P0 = P0*101325 #atm -> Pa

    VazaoVolumetrica =(sum(F)*9.437e-3)/ρ #mol/s -> m³/s
    v = VazaoVolumetrica/Ac #m³/s -> m/s
    vs = v*φ #velocidade superficial do gás

    ΔP_L = (150*μ*(1-φ**2))/(Dp**2*φ**3)*vs + 1.5*(vs**2*ρ*(1-φ))/(Dp*φ**3)
    P = P0 - L*ΔP_L
    P = P/101325 #Pa -> atm

    return P