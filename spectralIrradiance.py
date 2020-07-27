# WRITTEN BY MICHELE MARRO (michele.marro@unige.ch)

# THIS FILE CONTAINS THE CODE THAT REFERS TO
# BIRD ET AL 1986: SIMPLE SOLAR SPECTRAL MODEL
# DIRECT PART.

# CLOUD COVER CORRECTION FROM Jasmine et al, 1998
# to validate here

# FILE NEEDED: extra_irr, ozone_abs, unif_abs,
# water_abs, wl_list



from math import exp, pi, cos, sin
import matplotlib.pyplot as plt
import numpy as np
import csv

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


class SpectralIrradiance:
    
    # standard pressure at surface
    P0 = 1013.

    # current sub-division of external data
    # from 
    lat_grid = [89.95 - i/10 for i in range(1800)]
    lon_grid = [-179.95 + i/10 for i in range(3600)]

    def __init__(self, oZ = 0.28, preWater = 2., aers_opt_dep = 0.12, lat=0, lon=0, cloud_cover=False):
        """
        oZ : float
            ozone column in cm

        preWater : float 
            precipitable water in cm

        aers_opt_dep : float    
            aerosol optical depth at 0.5 micron
            Note: from 0.01 (clean atm) to 0.4
             -> mean value is 0.1-0.15
        
        lat : float
            latitude in degrees
            (needed for cloud cover)

        lon : float
            longitude in degrees
            (needed for cloud cover)
        """

        self.oZ = oZ
        self.preWater = preWater
        self.delta = aers_opt_dep

        # extra terrestrial irradiance
        self.E0 = extra_irr.data

        # ozone absoption coefficients
        self.k_o = ozone_abs.data

        # uniform mixed gases absoption coefficients
        self.k_g = unif_abs.data

        # water vapour absoption coefficients
        self.k_w = water_abs.data

        # wavelenghts used in this model
        self.wl = wl_list.data

        self.lat = lat
        self.lon = lon

        self.cloud_cover = cloud_cover

        if(self.cloud_cover):

            lat_index = min(range(len(self.lat_grid)), key=lambda i: abs(self.lat_grid[i] - self.lat))
            lon_index = min(range(len(self.lon_grid)), key=lambda i: abs(self.lon_grid[i] - self.lon))

            #if (len(self.E0), len(self.k_o), len(self.k_g), len(self.k_w)])) == 122):
            #    raise TypeError("Something gone wrong in external file dimension")
            
            # rows: lat
            # columuns: lon
            with open("cloud_cover/MODAL2_M_CLD_FR_2019-01-01_rgb_3600x1800_APR.CSV", mode='r') as csv_file:
                data = np.array([i for i in csv.reader(csv_file, delimiter=",",
                                            quoting=csv.QUOTE_NONNUMERIC)])

            self.CF = data[lat_index, lon_index]
        
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
        
        if(self.cloud_cover):
            cc_factor = [0.76 + 0.24*self.CF + 0.24*(1 - self.CF)*(i/0.49)**4 for i in self.wl]
        else:
            cc_factor = [1. for i in range(len(self.wl))] 
            self.CF = 1.

        # spectral power in W/m2/micron

        spectralIrradiance = [self.CF*cc*es_dist(day)*i*j*k*l*m*n \
        for cc, i, j, k, l, m, n in zip(cc_factor, self.E0, tau_r, tau_o, tau_g, tau_w, tau_a)]
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
        
        if(self.cloud_cover):
            cc_factor = [0.76 + 0.24*self.CF + 0.24*(1 - self.CF)*(i/0.49)**4 for i in self.wl]
        else:
            cc_factor = [1. for i in range(len(self.wl))]
            self.CF = 1.

        # spectral power in W/m2/micron

        spectralIrradiance = [self.CF*cc*es_dist(day)*i*j*k*l*m*n \
        for cc, i, j, k, l, m, n in zip(cc_factor, self.E0, tau_r, tau_o, tau_g, tau_w, tau_a)]

        #plt.plot(wl, tau_r, label="tau_r")
        #plt.plot(wl, tau_o, label="tau_o")
        #plt.plot(wl, tau_g, label="tau_g")
        #plt.plot(wl, tau_w, label="tau_w")

        plt.plot(self.wl,spectralIrradiance) 

        plt.title("ozone: {} cm - pre_water: {} cm - aerosol opt: {} - SZA: {} - P: {} - day : {}".format(self.oZ, self.preWater, self.delta, Z0, P, day))
        plt.show()

if __name__== "__main__":
    curr_irr = SpectralIrradiance(cloud_cover=True)
    curr_irr.get_irradiance(0, 1014, 1)
    curr_irr.plot_irradiance(0, 1014, 1)


