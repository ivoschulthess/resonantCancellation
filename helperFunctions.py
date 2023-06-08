# imports
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.constants import physical_constants
import matplotlib.pyplot as plt
import warnings
  
# suppress warnings
warnings.filterwarnings('ignore')


##############################
# CONSTANTS
##############################

# neutron de Broglie wavelength in [m]
lambda_n = 4.96e-10

# neutron mass in [kg]
m_n = physical_constants['neutron mass'][0]

# Planck constant in [J/Hz]
h = physical_constants['Planck constant'][0]

# neutron velocity in [m/s]
v_n = h / m_n / lambda_n

# gyromagnetic ratio of the neutron in [rad/s/T]
gamma_n = physical_constants['neutron gyromag. ratio'][0]


##############################
# FUNCTIONS
##############################

def weighted_mean(values, stds, axis=0):
    """
    Weighted mean and error of the mean.
    """
    average = np.average(values, weights=1/stds**2, axis=axis)
    std = 1 / np.sqrt(np.sum(1/stds**2, axis=axis))
        
    return np.array([average, std])

def ratio (A, B):
    """
    Ratio of two numbers/arrays and its error.
    """
    ratio = A[0] / B[0]
    ratioErr = np.sqrt((A[1]/B[0])**2 + (A[0]*B[1]/B[0]**2)**2)
    
    return np.array([ratio, ratioErr])

def sinFct(t, f, a, p, o):
    """
    Sinusoidal function.  
    """
    return abs(a) * np.sin(2*np.pi*f*t - p*np.pi/180) + o

def sinFrqFct (f, fp, f0, a, o):
    """
    Sinusoidal function, used to fit the Ramsey frequency scans. 
    """
    return abs(a) * np.sin(2*np.pi*fp*(f-f0)) + o

def sinPhaseFct (p, a, p0, o):
    """
    Sinusoidal function with a fixed period of 360, used to fit the Ramsey phase scans. 
    """
    return abs(a) * np.sin(np.pi/180*(p-p0)) + o

def butterworthFct(f, aMax, fCut, n=1):
    """
    higher-order low-pass filter function
    """
    return aMax / np.sqrt(1+(f/fCut)**(2*n))
    
def resonantCancelation(f, B, t_int, gamma=gamma_n): 
    """
    Resonant cancelation of the signal if an oscillating field is added in a Ramsey setup
    """
    return abs(gamma*B/np.pi/f * np.sin(np.pi*f*t_int))