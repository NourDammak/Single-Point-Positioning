# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 17:50:20 2021

@author: Nour Dammak
"""

import numpy as np
import math
from functions import *


nav_file = "almanac.yuma.week0131.405504.txt"
nav = read_yuma(nav_file)

#               year month day hour minute second
start_date = [  2021, 10,   10, 0,    0,      0]
stop_date =  [  2021, 10,   10, 24,    0,      0]


start_time = date2gpstime(start_date)
stop_time = date2gpstime(stop_date)

dt_5min = 5*60
dt_10min=10*60
dt_15min=15*60


nav1sat = nav[0,:]
nav2sat = nav[1,:]
nav3sat = nav[2,:]
""".
.
.
.
.
.
"""
nav32sat = nav[-1,:]



XYZ_ALL = np.zeros((0,4))


""" coordinates X Y Z for the sat 1 (PRN-01) evry 5 min"""
for t in range(start_time[1],stop_time[1],dt_5min):
    XYZ = satpos_alm(start_time[0],t,nav1sat)
    XYZ_T = np.hstack((t,XYZ))
    XYZ_ALL=np.append(XYZ_ALL,XYZ_T.reshape((1,4)), axis=0)
    # print t as this form "nbr of hours : number of minutes"
    # print(XYZ)
    print(t//3600,":",t%3600//60 ,XYZ)




