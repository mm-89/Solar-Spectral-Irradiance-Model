# WRITTEN BY MICHELE MARRO (michele.marro@unige.ch)

# THIS FILE CONTAINS THE CODE THAT REFERS TO
# BIRD ET AL 1986: SIMPLE SOLAR SPECTRAL MODEL
# DIRECT PART.

# FILE NEEDED: extra_irr, ozone_abs, unif_abs,
# water_abs, wl_list



from math import exp, pi, cos, sin
import matplotlib.pyplot as plt
import numpy as np

try:

    import extra_irr
    import ozone_abs
    import unif_abs
    import water_abs
    import wl_list

except:
    raise TypeError("PROBLEM TO IMPORT OR FIND EXTERNAL DATA")

# ----------------
# FUNCTIONS

def m_p(Z0, P0, P):
    return m(Z0) * P/P0

def m_0(Z0):
    return 35./( 1224. * cos(Z0)**2  + 1 )**0.5

def m(Z0):
    return 1/( cos(Z0) + 0.15*(93.885 - Z0)**(-1.253) )

def day_angle(day):
    """
    day from 1 to 365
    """
    return 2*pi*( day - 1 )/365

def es_dist(day):
    return 1.00011 + \
            0.034221*cos(day_angle(day)) + \
            0.00128*sin(day_angle(day)) + \
            0.000719*cos(2*day_angle(day)) + \
            0.000077*sin(2*day_angle(day))


class SpectralIrradiance:
    
    # standard pressure at surface
    P0 = 1013.

    def __init__(self, oZ = 0.28, preWater = 2., aers_opt_dep = 0.12):
        """
        oZ : float
            ozone column in cm

        preWater : float 
            precipitable water in cm

        aers_opt_dep : float    
            aerosol optical depth at 0.5 micron
            Note: from 0.01 (clean atm) to 0.4
             -> mean value is 0.1-0.15
        """
        self.oZ = oZ
        self.preWater = preWater
        self.delta = aers_opt_dep

        self.E0 = extra_irr.data
        self.k_o = ozone_abs.data
        self.k_g = unif_abs.data
        self.k_w = water_abs.data

        self.wl = wl_list.data

        #if(np.any([len(self.E0), len(self.k_o), len(self.k_g), len(self.k_w)]) != 122):
        #    raise TypeError("Something gone wrong in external file dimension")
        
        print("***")
        print("Simple solar spectral irradiance model")
        print("***")
        print("Values are: ")
        print(" ozone column in cm: {}".format(self.oZ))
        print(" precipitable water in cm: {}".format(self.preWater))
        print(" aerosol optical depth at 0.5 micron: {}".format(self.delta))
        print("***")

    def get_irradiance(self, Z0, P, day):

        Z0 = Z0*pi/180.
                
        tau_r = [exp( -m_p(Z0, self.P0, P) / ( (i**4) * (115.6406 - 1.335/i**2)) ) for i in self.wl]
        tau_o = [exp( -i * self.oZ * m_0(Z0) ) for i in self.k_o]
        tau_g = [exp(-1.41 * i * m_p(Z0, self.P0, P) / (1 + 118.3 * i * m_p(Z0, self.P0, P))**(0.45) ) for i in self.k_g]
        tau_w = [exp( -0.2385 * i * self.preWater * m(Z0)/(1 + 20.07 * i* self.preWater * m(Z0))**(0.45) ) for i in self.k_w]
        tau_a = []
        for i in self.wl:
            if(i<0.5):
                tau_a.append( exp(-self.delta*(i/0.5)**(-1.0274)*m(Z0)) )
            else:
                tau_a.append( exp(-self.delta*(i/0.5)**(-1.2060)*m(Z0)) )

        # spectral power in W/m2/micron

        spectralIrradiance = [es_dist(day)*i*j*k*l*m*n for i, j, k, l, m, n in zip(self.E0, tau_r, tau_o, tau_g, tau_w, tau_a)]
        return spectralIrradiance


    def plot_irradiance(self, Z0, P, day):

        Z0 = Z0*pi/180.
                
        tau_r = [exp( -m_p(Z0, self.P0, P) / ( (i**4) * (115.6406 - 1.335/i**2)) ) for i in self.wl]
        tau_o = [exp( -i * self.oZ * m_0(Z0) ) for i in self.k_o]
        tau_g = [exp(-1.41 * i * m_p(Z0, self.P0, P) / (1 + 118.3 * i * m_p(Z0, self.P0, P))**(0.45) ) for i in self.k_g]
        tau_w = [exp( -0.2385 * i * self.preWater * m(Z0)/(1 + 20.07 * i* self.preWater * m(Z0))**(0.45) ) for i in self.k_w]
        tau_a = []
        for i in self.wl:
            if(i<0.5):
                tau_a.append( exp(-self.delta*(i/0.5)**(-1.0274)*m(Z0)) )
            else:
                tau_a.append( exp(-self.delta*(i/0.5)**(-1.2060)*m(Z0)) )

        # spectral power in W/m2/micron

        spectralIrradiance = [es_dist(day)*i*j*k*l*m*n for i, j, k, l, m, n in zip(self.E0, tau_r, tau_o, tau_g, tau_w, tau_a)]

#plt.plot(wl, tau_r, label="tau_r")
#plt.plot(wl, tau_o, label="tau_o")
#plt.plot(wl, tau_g, label="tau_g")
#plt.plot(wl, tau_w, label="tau_w")

        plt.plot(self.wl,spectralIrradiance) 

        plt.title("ozone: {} cm - pre_water: {} cm - aerosol opt: {} - SZA: {} - P: {} - day : {}".format(self.oZ, self.preWater, self.delta, Z0, P, day))
        plt.show()




