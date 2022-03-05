# -*- coding: utf-8 -*-
"""
Created on Tue May  5 13:42:38 2020

@author: Maciek
"""
import numpy as np
from datetime import date
import math as mat

C=299792458
# MU=3.986004415e14
MU = 3.986005e14
OMGE=7.2921151467e-5
# MU        = 3.986005e14       # WGS84 Earth's gravitational constant for GPS user [m^3/s^2]
# OMGE    = 7.2921151467e-5

def satpos(week, tow, nav):
    # czas tygodnia GPS, a nie doby
    toe =  nav[:, 17]
    t_full = week * 604800 + tow
    toe_full = nav[:,27] * 604800 + toe
    tk = t_full - toe_full
    i=np.argmin(abs(tk))
    tk=tk[i]
    toe = toe[i]
    nav=nav[i,:]
    Asqrt = nav[16]
    dn  =  nav[11]
    M0  =  nav[12]
    e   =  nav[14] 
    omg =  nav[23]
    cus =  nav[15]
    cuc =  nav[13]
    cis =  nav[20]
    cic =  nav[18]
    crs = nav[10]
    crc =  nav[22]
    i0  =  nav[21]
    omgdot = nav[24]
    di_dt = nav[25]
    L0 = nav[19]
    af0 = nav[6]
    af1 = nav[7]
    af2 = nav[8]
    
    A = Asqrt**2                #wielka półos orbity                  #mimosród
    n = mat.sqrt(MU/A**3) + dn #ruch redni (wyznaczony i poprawiony)
    Mk = M0+n*tk              #anomalia srednia
    Ek = Mk
    while 1:
        Eprev = Ek
        Ek = Mk + e*np.sin(Eprev)
        if (abs(Ek-Eprev)<1e-12): break       #anomalia mimosrodowa
    
    # Ek = Mk + (e * mat.sin(Mk))/(1- mat.sin(Mk + e) + mat.sin(Mk))
    # while 1:
    #     Eprev = Ek
    #     Ek = Eprev - (Eprev - e*mat.sin(Eprev) - Mk)/(1-e*mat.cos(Eprev))
    #     print(Ek)
    #     if (abs(Ek-Eprev)<1e-15): break 
    
    vk = mat.atan2(mat.sqrt(1-e**2)*mat.sin(Ek),mat.cos(Ek)-e)  #anomalia prawdziwa
    # argument szerokosci
    fik = vk + omg
    
    #poprawki do argumentu szerokosci, promienia orbity i inklinacji
    duk = cus * mat.sin(2*fik) + cuc * mat.cos(2*fik)
    drk = crs * mat.sin(2*fik) + crc * mat.cos(2*fik)
    dik = cis * mat.sin(2*fik) + cic * mat.cos(2*fik)
    
    # r1 = np.array(([uk],[A*(1-e*mat.cos(Ek))],[i0+di_dt*tk])) #
    # r2 = np.dot(np.array(([cus,cuc],[crs,crc],[cis,cic])),np.array(([mat.sin(2*uk)],[mat.cos(2*uk)]))) #poprawki perturbacyjne do argumentu szerokosci
    # r = r1+r2                     #poprawienie argumentu szerokosci, wektora wodzącego satelity i nachylenia orbity
    # poprawione argument szerokosci, promień orbity i inklinacja
    uk = fik + duk
    rk = A * (1 - e*mat.cos(Ek)) + drk
    ik = i0 + di_dt * tk + dik
    
    OMGk = L0 + (omgdot - OMGE)*tk - OMGE*toe #długosc wezła wstępującego
    xk = rk * mat.cos(uk)
    yk = rk * mat.sin(uk)

    Xk = xk*mat.cos(OMGk) - yk*mat.cos(ik) * mat.sin(OMGk)
    Yk = xk*mat.sin(OMGk) + yk*mat.cos(ik) * mat.cos(OMGk)
    Zk = yk*mat.sin(ik)
    
    XYZ = np.array([Xk,Yk,Zk])
    
    if nav[0]>1900:
        y = 0
    elif nav[0]<70:
        y=2000
    else:
        y=1900
    # obliczenie dt
    # week_alfa,t_alfa = date2tow([int(nav[0]+y), int(nav[1]),int(nav[2]), int(nav[3]), int(nav[4]), int(nav[5])], rollover = False)
    # talfa_full = week_alfa * 604800 + t_alfa
    # tc = t_full - talfa_full
    # dt = np.dot(nav[6:9],[1,tc,tc**2]) -2*mat.sqrt(MU*A)*e*mat.sin(Ek)/C**2
    dt = af0 + af1 * tk + af2 * tk**2 - (2*mat.sqrt(MU)/C**2)*e*mat.sqrt(A)*mat.sin(Ek)
    return XYZ, dt

a = 6378137
e2 = 0.00669438002290

def hirvonen(X,Y,Z):
    r = (X**2 + Y**2)**0.5
    B = mat.atan(Z/(r*(1-e2)))
    
    while 1:
        N = Np(B)
        H = r/np.cos(B) - N
        Bst = B
        B = mat.atan(Z/(r*(1-(e2*(N/(N+H))))))    
        if abs(Bst-B)<(0.00001/206265):
            break
    L = mat.atan(Y/X)
    N = Np(B)
    H = r/np.cos(B) - N
    return B, L, H

def Np(B):
    N = a/(1-e2*(np.sin(B)**2))**0.5
    return N


def tropo(el,h):
    
    p0 = 1013.25
    t0=291.5
    Rh0=0.5
    p=p0*(1-(0.0000226*h))**5.225
    t=t0-0.0065*h
    Rh=Rh0*mat.exp(-0.0006396*h)
    e=6.11*Rh*10**(7.5*(t-273.15)/(t-35.85))
    Tdz=0.002277*p
    twz=0.002277*(1255/(t+0.05))*e   
    trop=1/(np.sin(el))*(Tdz+twz)
    return trop

    
def date2gpstime(date1):
    
    dif_of_days = date.toordinal(date(date1[0], date1[1],date1[2])) - date.toordinal(date(2019,4,7))
    
    weeks = dif_of_days//7
    
    day_of_week = dif_of_days-weeks*7
    
    #tow = time of the week
    
    tow = day_of_week*86400 + date1[3]*3600 + date1[4]*60 + date1[5]
     
    return [weeks, tow]


