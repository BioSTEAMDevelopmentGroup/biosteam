""" A module of tables and graphs for constants
"""
from math import log, pi

__all__ = ('GTable', 'HNATable', 'ceil_half_step',
           'compute_vessel_weight_and_wall_thickness',
           'compute_Stokes_law_York_Demister_K_value')

def GTable(DRho, Hlr):
    """
    Return the allowable downflow (baffle liquid load) in gph/ft2, usually used for vertical vessel.
    
    Parameters
    ----------
    DRho : float
        Density difference between light liquid and vapor [lb/ft^3 ?]
    Hlr : float
        Height of liquid level above the interphase of light liquid and heavy liquid [ft]
    
    """

    A = {}
    B = {}
    C = {}
    D = {}

    # TODO: Add errors when outside the range!
    if Hlr > 30.0:
        Hlr = 30

    if Hlr < 18.0:
        Hlr = 18

    if DRho > 50.0:
        DRho = 50

    if DRho < 10.0:
        DRho = 10
    
    Hlr = round(Hlr, 0)

    A[18] = -9000.0
    B[18] = 1275.4
    C[18] = -31.3571
    D[18] = 0.255556

    A[19] = -4690.0
    B[19] = 900.117
    C[19] = -20.5252
    D[19] = 0.157254

    A[20] = -9980.0
    B[20] = 1367.91
    C[20] = -33.0163
    D[20] = 0.26476

    A[21] = -8120.0
    B[21] = 1147.46
    C[21] = -25.3
    D[21] = 0.184444

    A[22] = -16800.0
    B[22] = 1964.89
    C[22] = -48.8627
    D[22] = 0.399498

    A[23] = -7900.0
    B[23] = 1255.35
    C[23] = -29.9142
    D[23] = 0.235632

    A[24] = -11200.0
    B[24] = 1561.48
    C[24] = -38.7335
    D[24] = 0.318511

    A[25] = -11100.0
    B[25] = 1554.66
    C[25] = -38.0313
    D[25] = 0.308026

    A[26] = -7410.0
    B[26] = 1274.0
    C[26] = -30.8013
    D[26] = 0.246585

    A[27] = -12700.0
    B[27] = 1709.78
    C[27] = -42.1048
    D[27] = 0.342222

    A[28] = -10200.0
    B[28] = 1507.78
    C[28] = -36.422
    D[28] = 0.291221

    A[29] = -10700.0
    B[29] = 1553.51
    C[29] = -37.5721
    D[29] = 0.300279

    A[30] = -9830.0
    B[30] = 1513.11
    C[30] = -37.1907
    D[30] = 0.30379

    G = A[Hlr] + (B[Hlr] * DRho) + (C[Hlr] * DRho ** 2) + (D[Hlr] * DRho ** 3)

    return round(G, 2)


def HNATable(Type, X):
    """
    Table for cylindrical height and area conversions.
        
    Parameters
    ----------
    Type: int
        1 if given H/D and find A/At, 2 if given A/At and find H/D.
    X: float
        H/D or A/At.
        
    """
    # Type = 1 is where H/D is known, find A/At, Type = 2 is where A/At is known, find H/D
    if (Type == 1):
        a = -0.0000475593
        b = 3.924091
        c = 0.174875
        d = -6.358805
        e = 5.668973
        f = 4.018448
        g = -4.916411
        h = -1.801705
        i = -0.145348
        Y = (a + c * X + e * X ** 2 + g * X ** 3 + i * X ** 4) / \
            (1.0 + b * X + d * X ** 2 + f * X ** 3 + h * X ** 4)
    else:
        a = 0.00153756
        b = 26.787101
        c = 3.299201
        d = -22.923932
        e = 24.353518
        f = -14.844824
        g = -36.999376
        h = 10.529572
        i = 9.892851
        Y = (a + c * X + e * X ** 2 + g * X ** 3 + i * X ** 4) / \
            (1.0 + b * X + d * X ** 2 + f * X ** 3 + h * X ** 4)

    return Y


def compute_vessel_weight_and_wall_thickness(P, D, L, rho_M, Je=0.85):
    """Return vessel weight and wall thickness.
    
    Parameters
    ----------
    P : float 
     Pressure [psia].
    D : float
        Diameter [ft].
    L: float
        Vessel length [ft].
    rho_M: float
        Density of Material [lb/ft^3].
    Je: float
        Joint efficiency (1.0 for X-Rayed joints, 0.85 for thin carbon steel),
    
    """
    S = 15000.0     # Vessel material stress value (assume carbon-steel)
    Ca = 1.0/8.0    # Corrosion Allowance in inches
    
    P_gauge = abs(P - 14.7)
    P1 = P_gauge + 30.0
    P2 = 1.1 * P_gauge
    if P1 > P2:
        PT = P1
    else:
        PT = P2

    # Calculate the wall thickness and surface area
    # Shell
    SWT = (PT * D*12.0) / (2.0 * S * Je - 1.2 * PT) + Ca
    SSA = pi * D * L
    if D < 15.0 and PT > (100 - 14.7):
        # Elliptical Heads
        HWT = (PT * D*12.0) / (2.0 * S * Je - 0.2 * PT) + Ca
        HSA = 1.09 * D ** 2
    elif D > 15.0:
        # Hemispherical Heads
        HWT = (PT * D*12.0) / (4.0 * S * Je - 0.4 * PT) + Ca
        HSA = 1.571 * D ** 2
    else:
        # Dished Heads
        HWT = 0.885 * (PT * D*12.0) / (S * Je - 0.1 * PT) + Ca
        HSA = 0.842 * D ** 2

    # Approximate the vessel wall thickness, whichever is larger
    if SWT > HWT:
        ts = SWT
    else:
        ts = HWT

    # Minimum thickness for vessel rigidity may be larger
    if D < 4:
        ts_min = 1/4
    elif D < 6:
        ts_min = 5/16
    elif D < 8:
        ts_min = 3/8
    elif D < 10:
        ts_min = 7/16
    elif D < 12:
        ts_min = 1/2
    else:
        ts_min = ts
    if ts < ts_min:
        ts = ts_min
    VW = rho_M * ts/12 * (SSA + 2.0 * HSA)  # in lb
    VW = round(VW, 2)
    return VW, ts


def compute_low_liq_level_height(Type, P, D):
    """Table to obtain Hlll for two-phase separators.
    
    Parameters
    ----------
    Type : int
        1 for vertical, 2 for horizontal.    
    P : float
        Pressure [psia].
    D: float 
        Diameter [ft].
    
    """
    if Type == 1:
        Hlll = 0.5
        if P < 300:
            Hlll = 1.25

    elif Type == 2:
        if D <= 4.0:
            Hlll = 9.0/12.0
        elif D > 4.0 and D <= 7.0:
            Hlll = 10.0/12.0
        elif D > 7.0 and D <= 9.0:
            Hlll = 11.0/12.0
        elif D > 9.0 and D <= 11.0:
            Hlll = 1.0
        elif D > 11.0 and D <= 15.0:
            Hlll = 13.0/12.0
        else:
            Hlll = 15.0/12.0

    return Hlll  # in ft


def compute_Stokes_law_York_Demister_K_value(P):
    """
    Return K-constant in Stoke's Law using the York-Demister equation.
    
    Parameters
    ----------
    P : float
        Pressure [psia].
        
    """
    if P >= 0 and P <= 15.0:
        K = 0.1821+(0.0029*P)+(0.046*log(P))
    elif P > 15.0 and P <= 40.0:
        K = 0.35
    elif P > 40.0 and P <= 5500.0:
        K = 0.43 - 0.023*log(P)
    else:
        raise ValueError(f'invalid Pressure {P} psia')
    return K


def ceil_half_step(value):
    """Return value to the next/highest 0.5 units"""
    intval = round(value)
    if value > intval:
        return intval + 0.5
    elif value == intval:
        return value
    else:
        return intval - 0.5
