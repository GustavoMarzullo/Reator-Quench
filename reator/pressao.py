import numpy as np

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