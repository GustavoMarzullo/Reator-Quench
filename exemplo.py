from Reator.reator import Reagente, Leito, Reator

reagente = Reagente(
    FN2=0.21825*6500,
    FH2=0.65475*6500,
    FNH3=0.05*6500,
    P=180,
    T=400)

L1 = Leito(
    L = 5,
    phi=0.4,
    Ac=7,
    Dp=2)

L2 = Leito(
    L = 15,
    T = 380,
    phi=0.4,
    Ac=7,
    Dp=2)

L3 = Leito(
    L = 15,
    T = 380,
    phi=0.4,
    Ac=7,
    Dp=2)

reator = Reator(
    Leitos=[L1, L2, L3],
    Reagente=reagente,
    Y = [0.5, 0.3, 0.2])

resultado = reator.solve()

resultado.plot()
