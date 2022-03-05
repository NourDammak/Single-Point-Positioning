# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 16:50:58 2022

@author: HP
"""

import math as m
import numpy as np

def tropo(el,h):
    
    p0 = 1013.25
    t0=291.5
    Rh0=0.5
    p=p0*(1-(0.0000226*h))**5.225
    t=t0-0.0065*h
    Rh=Rh0*m.exp(-0.0006396)
    e=6.11*Rh*10**(7.5*(t-273.15)/(t-35.85))
    Tdz=0.0027*p
    twz=0.0027*(1255/(t+0.05))*e   
    trop=1/(np.sin(el))*(Tdz+twz)
    return print(trop)

tropo(50,30)