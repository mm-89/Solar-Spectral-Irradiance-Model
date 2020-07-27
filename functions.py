from math import pi, cos, sin

# ----------------
# FUNCTIONS

def m_p(Z0, P0, P):
    """
    Pressure-corrected 
    relative air mass.
    """
    return m(Z0) * P/P0

def m_0(Z0):
    """
    Relative air-mass
    for ozone.
    
    From: Paltridge and Platt, 1976
    """
    return 35./( 1224. * cos(Z0)**2  + 1 )**0.5

def m(Z0):
    """
    Relative air-mass.

    From: Kasten, 1966
    """
    return 1/( cos(Z0) + 0.15*(93.885 - Z0)**(-1.253) )

def day_angle(day):
    """
    Function that converts 
    day in a year in radiant
    -----------------
    Parameter:
    day : int
        from 1 to 365
    """
    return 2*pi*( day - 1 )/365

def es_dist(day):
    """
    Correction factor for Eath-Sun distance.

    From: Spencer, 1971
    """
    return 1.00011 + \
            0.034221*cos(day_angle(day)) + \
            0.00128*sin(day_angle(day)) + \
            0.000719*cos(2*day_angle(day)) + \
            0.000077*sin(2*day_angle(day))
