from Reator.reator import Reagente, Leito, Reator

Ac = 0.95
F= 6477 #mol/s
reagente = Reagente(
    FN2=0.2173*F,   
    FH2=0.6608*F,
    FNH3=0.0380*F,
    P=286,
    T=427)

L1 = Leito(
    L = 4.07/Ac,
    Ac = 0.95)

L2 = Leito(
    L = 4.07/Ac,
    Ac = 0.95,
    T = 327)

reator = Reator(
    Leitos=[L1, L2],
    Reagente=reagente,
    Y = [0.5, 0.5])

resultado = reator.solve()

resultado.plot()