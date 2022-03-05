# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 23:48:08 2022

@author: HP
"""
import numpy as np
from pylab import *
import math


#klobucgar Function

    
alpha = [1.2107e-08, -7.4506e-09, -1.1921e-07, 5.9605e-08]

beta = [9.8304e+04, -8.1920e+04, -1.9661e+05, 4.5875e+05]
    

tpgs = 0 #epch of time

#el=30
#az=180
#fi=50
#lam=20


def klobuchar(el,az,fi,lam,alfa,beta,tgps):
    
    alpha = [1.2107e-08, -7.4506e-09, -1.1921e-07, 5.9605e-08]

    beta = [9.8304e+04, -8.1920e+04, -1.9661e+05, 4.5875e+05]
    
    els=el/180
    azs=az/180
    fis=fi/180
    lams=lam/180
    
    
    # calculate Earch center
    # 1. Earth centered angla of IPP
    
    psi=0.0137/(els+0.11)-0.022
    print(psi)
    
    # 2.latitudeof IPP
    
    fi_ipp=fis+ psi * np.cos(np.deg2rad(az))
    
    if fi_ipp > 0.416:
        fi_ipp = 0.416
    elif fi_ipp<-0.416:
        fi_ipp = -0.416
        
        
    # 3. Longitude of IPP
    
    lam_ipp = lams + (psi*np.sin(np.deg2rad(az))/np.cos(fi_ipp*np.pi))
    print(lam_ipp)
    
    # 4 Geomagnetic latitude of IPP
    fi_mag_ipp = fi_ipp + 0.064 * np.cos((lam_ipp - 1.617)*np.pi)
    
    print(fi_mag_ipp)
    
    # 5 Localtime
    
    tsow = 43200 * lam_ipp + tgps
    print(tsow)
    
    tsod = math.fmod(tsow,86400)
    
    if tsod>68400:
        tsod = tsod-86400
    elif tsod<0:
        tsod = tsod + 86400
        
    print(tsod)
    
    #6 Amplitude of ionospheric delay
    
    aion = alpha[0] + alpha[1] * fi_mag_ipp + alpha[2]*fi_mag_ipp**2 + alpha[3]*fi_mag_ipp**3
    if aion<0:
        aion=0
    print(aion)
    
    #7 Period ofionospheric delay
    
    pion = beta[0]+ beta[1]*fi_mag_ipp + beta[2]*fi_mag_ipp**2 + beta[3]*fi_mag_ipp**3
    if pion<72000:
        pion = 72000
    print(pion)
    
    #8 Phase of the iono delay
    
    fi_ion=2*np.pi*(tsod-50400)/pion
    print(fi_ion)
    
    
    #9 slant factor
    
    mf=1+16*(0.53 - els)**3
    
    print(mf)
    
    #Actual ionospheric delay value
    
    C = 299792458.0
    
    if abs(fi_ion)<np.pi/2:
        dI = C*mf*(5*10**-9 + aion*(1-(fi_ion**2)/2 + (fi_ion**4)/24))
    elif abs(fi_ion)>np.pi/2:
        dI=C*mf*5*10**-9
    
    return (dI)
    
    
    
    
    
    
    


