# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:15:49 2021

@author: Maciek
SPP z odległoscią geometryczną w pętli
"""
from readrnx import readrnxobs, readrnxnav, date2tow
from functions import satpos, hirvonen, tropo
import numpy as np
import math as m 
from klobuchar import  klobuchar
# path to the broadcast navigation data - navigation RINEX file
nav_file = "BRDC00WRD_R_20213180000_01D_GN.rnx"
# path to the observation data - observation RINEX file
obs_file = "KRA100POL_R_20213180000_01D_30S_MO.rnx"
  
# define time of observation file - should be the time the same that in the observation file
# for first epoch only:
time_start = [2021, 11, 14,  0,  0,  0]
time_end = [2021, 11, 14,  0,  0,  0]

#convert the epoch to number of second (we need week and tr)
week, tr, kl = date2tow(time_start)

# reading the observation data from the file
obs,iobs = readrnxobs(obs_file, time_start, time_end, 'G')
# reading the navigation data from the file:
nav, inav = readrnxnav(nav_file)
#%%
# defining the approximate receiver coordinates: X,Y,Z
Xtrue = np.array([3856937.5695,  1397749.6854,  4867719.1135])
Xr0 = np.array([3850000.0000,  1390000.0000,  4860000.0000])
# some settings:
el_mask = 10 # elevation mask/cut off in degrees
# for the calculations only of the first epoch in the file, t is equal to 0 (because the data are from Sunday),
# week equal to 2184. If uou want to make calculations for different epochs, you can use the function date2tow(data), which can
# be found in the readrnx file, using the same as in the project 1

e2= 0.00669438002290
a=6378137


n_sat_all=[]
GDOP_0=[]
n_sat=0
A11=np.zeros((0,4))


#taking only the observations from the specified epoch
i = iobs[:,2]==tr
pobs=obs[i,:]
iobst = iobs[i,:]


#klobucgar Function : in klobichar.py file
alpha = [1.2107e-08, -7.4506e-09, -1.1921e-07, 5.9605e-08]
beta = [9.8304e+04, -8.1920e+04, -1.9661e+05, 4.5875e+05]
# observation from one epoch - Pseudorange obs
#PRNs of satellites from the epoch
sats = iobst[:,0]
dtr = 0
tau=0.07
# create a vector of ro (geometric distance) for all satellites, in the first epoch ro can be tau*C, where tau is approximate transmission time, equal to 0.07s
# in the first epoch tau should be 0.07 and ro should be tau*C

Xs_tab=np.zeros((int(nav[-1,0])+1, len(np.arange(5)),5))*np.nan
Vto=np.ones((10))*0.07
Xr=Xr0
for i in range (1):
    # loop for all the satellites in the epoch, good ideas is to get also index, using enumerate
    A=np.zeros((0,4)) 
    Vro=np.zeros((0,10))
    Trop_Values=np.zeros((0,8))
    GDOP_0=[]
    PDOP_0=[]
    TDOP_0=[]
    PDOPneu_0=[]
    HDOP_0=[]
    VDOP_0=[]
    el_values=[]
    Az_values=[]
    Vpobs=np.zeros((8,0))
    Iono_Values=np.zeros((8,0))
    Pcalc_Values=np.zeros((8,0))
    Vector_L=[]
    lsat=0
    ind=[]
    Xr_all=np.zeros((0,3))
    
    
    for isat, sat in enumerate(sats):
        # calculate transmission epoch: ts = tr - tau + dtr, dtr is the receiver clock error, in the first iteration, it is zero
        ts = tr - Vto[isat] + dtr
        # taking the ephemerides for one satellite - it should be an array of nav data, use (inav==sat)
        nav1sat = nav[(inav==sat),:]
        
        # calculating the coordinates of satellites, for the transmission epoch, use the satpos function, see the input values in the code
        # the outputs should be the XYZ and dts: inputs: satpos(week, ts, nav1sat)
        Xs, dts = satpos(week,ts,nav1sat)
        # rotation of the satellite coordinates due to rotation of the earth, use the rotation around Z-axis
        w=0.000072921151467
        C=299792458
        
        Xs_Mat = np.array([[np.cos(w*Vto[isat]), np.sin(w*Vto[isat]), 0],[ -np.sin(w*Vto[isat]),np.cos(w*Vto[isat]), 0 ],[ 0,0,1]])
    
        Xs_rot=Xs_Mat.dot(Xs)
        # satellite-receiver vector
        Xsr = Xs_rot - Xr
              
        # the geometric distance between the satellite and the receiver
       
        prs=m.sqrt((Xsr[0])**2+ (Xsr[1])**2+(Xsr[2])**2 )
           
        Vto[isat]=prs/C
        ro = Vto[isat]*C
        
        Vro=np.append(Vro,ro, axis=None)
        
        # elevation and azimuth calculations
        B,L,H = hirvonen(Xtrue[0],Xtrue[1],Xtrue[2])
        phi=B
        lamda=L
        h=227.5428
        
        
        r=np.linalg.norm(Xs-Xr)
        nv=np.array([[-m.sin(phi)*m.cos(lamda),-m.sin(phi)*m.sin(lamda), m.cos(phi)]])
        ev=np.array([[-m.sin(lamda), m.cos(lamda),0]])
        uv=np.array([[m.cos(phi)*m.cos(lamda),m.cos(phi)*m.sin(lamda),m.sin(phi)]])
        R=np.hstack([nv.T, ev.T, uv.T])


        Xrneu=np.dot (R.T,Xsr.T)

        #Calculation of the azimuth(Az) and elevation(el) angles
        n=Xrneu[0]
        e=Xrneu[1]
        u=Xrneu[2]
        Az=np.rad2deg(m.atan2(e,n))
        Az_values=np.append(Az_values,Az, axis=None)
        el=np.rad2deg(m.atan2(u,m.sqrt(e**2+n**2)))
        el_values=np.append(el_values,el, axis=None)

        # if elevation above elevation mask, do the following calculations
        if el > el_mask:
            A1row=np.array([[-(Xs_rot[0]-Xr[0])/r, -(Xs_rot[1]-Xr[1])/r, -(Xs_rot[2]-Xr[2])/r,1 ]])
            A=np.append(A,A1row, axis=0)
            lsat+=1
            Xs_tab[int(sat),0,0]=el
            Xs_tab[int(sat),0,1]=Az
            Xs_tab[int(sat),0,2]=Xs[0]
            Xs_tab[int(sat),0,3]=Xs[1]
            Xs_tab[int(sat),0,4]=Xs[2]
        
        # calculate the atmospehric corrections: tropo and iono
        #Tropo
        if el > el_mask:
            elrad=m.radians(el) 
            Trop=tropo(elrad,h)
            Trop_Values=np.append(Trop_Values,Trop, axis=None)
        #Iono
            Iono=klobuchar(el, Az, phi, lamda, alpha, beta, tr)
            Iono_Values=np.append(Iono_Values,Iono, axis=None)
            
        # compute calculated Pseudorange: Pcalc = ro + C*dtr-C*dts + tropo + iono
        # in testing, ommit the atmosphere, and assume that tropo = 0 and iono = 0

            Pcalc = ro + C*dtr-C*dts + Iono + Trop
            Pcalc_Values=np.append(Pcalc_Values,Pcalc,axis=None)
            # compute vector L: L = Pcalc - Pobs
            L1 = Pcalc - obs[isat]
            Vector_L=np.append(Vector_L,L1,axis=None) 
            ind.append(isat)
            
            
            
        # Computation of number of sattelite and the  for this epoch      
        if el > el_mask:
            A1row=np.array([[-(Xs[0]-Xr[0])/r, -(Xs[1]-Xr[1])/r, -(Xs[2]-Xr[2])/r,1 ]])
            A11=np.append(A11,A1row, axis=0)
            n_sat+=1
    #DOP values
    Q = np.linalg.inv(np.dot(A11.T,A11))
    Qxyz=np.array([[Q[0,0],Q[0,1],Q[0,2]],[Q[1,0],Q[1,1],Q[1,2]],[Q[2,0],Q[2,1],Q[3,2]]])
    Qneu=np.dot(R.T,Qxyz,R)
    GDOP=np.sqrt(Q[0,0] + Q[1,1] + Q[2,2] + Q[3,3])
    PDOP=np.sqrt(Q[0,0] + Q[1,1] + Q[2,2])
    TDOP=np.sqrt(Q[3,3])
    PDOPneu=np.sqrt(Qxyz[0,0] + Qxyz[1,1] + Qxyz[2,2])
    HDOP=np.sqrt(Qxyz[0,0] + Qxyz[1,1])
    VDOP=np.sqrt(Qxyz[2,2])
    
    #to verify that PDOPneu equal to sqrt((HDOP**2+VDOP**2)))
    #a=PDOPneu-(np.sqrt((HDOP**2+VDOP**2)))
    #print("a=",a)
    
    GDOP_0.append(GDOP)
    PDOP_0.append(PDOP)
    TDOP_0.append(TDOP)
    PDOPneu_0.append(PDOP)
    HDOP_0.append(HDOP)
    VDOP_0.append(VDOP)
    n_sat_all.append(n_sat)
    a1=Pcalc_Values-pobs[ind,0]          
    # solve the equation system using least-squares method (after all satellites processed)
    # use function: np.dot, np.linalg.inv (A och FF)
    a2=np.dot(A.T,a1)
    a3=np.dot(A.T,A)
    a4=-np.linalg.inv(a3)
    x=np.dot(a4,a2)
    # modify the X, Y and Z using the results of the least square method
    Xr = Xr + x[0:-1]
    Xr_all=np.append(Xr_all,Xr, axis=None)
	# modify the receiver clock error, remember to divise the resulting cdt by speed of light (C) before adding to dt_r
    dtr = np.array(dtr+x[-1]/C) 

    # if all iterations are done, save the epoch, X, Y, Z and dt_r for some whole day analysis



#coordinates of the sattelite
Xs2, dts = satpos(week,ts,nav1sat)


#coordinates of receiver without corrections
N= a/(m.sqrt(1-(e2*m.sin(np.deg2rad(phi)**2))))

Xk = (N+h)*(m.cos(phi))*(m.cos(lamda))
Yk = (N+h)*m.cos(phi)*(m.sin(lamda))
Zk=(N*(1-e2)+h)*m.sin(phi)

#Xr : coordinates of the receiver without correction
XYZ_no_correction = np.array([Xk,Yk,Zk])





#
