#! /usr/bin/env python


import numpy as np
from ppclass import pp
import ppplot
import ppcompute
import matplotlib.pyplot as plt
import matplotlib.mlab as lab
import matplotlib.cm as cm
from scipy import signal 
from sympy.mpmath import *

vv="W"
tt=600
xx=75
xm="0,1e10"
zz=161
ww="W"
tp="T"
#fi="Venus4070_CKJan15H_0_minuit_0/wrfout_d01_9999-01-01_22:13:20"
#fi="Venus4070_CKJan15H_0_midi_100tes/wrfout_d01_9999-01-01_22:13:20"
fi="Venus4070_CKMai14Z_0_minuit/wrfout_d01_9999-01-02_03:46:40"
#fi="Venus4070_1d_0_minuit/wrfout_d01_9999-01-01_11:06:40"
fi="wrfout_d01_9999-01-03_07:33:20"
pl = ppplot.plot2d()
p1 = ppplot.plot1d()

p, x,y,z,t = pp(file=fi,var="PTOT",t=tt,x=xm,y=xm).getfd()
np.savetxt('pmoy.txt',p)
x1=85
x2=195
Wxy = pp(file=fi,var=ww,t=tt).getf()[zz]
Wm = pp(file=fi,var=ww,t=tt,x=xm,y=xm).getf()[zz]
Tz= pp(file=fi,var=tp,t=tt,y=xx,x=xx).getf()
Tm=pp(file=fi,var=tp,t=tt,x=xm,y=xm).getf()

CST=np.loadtxt('filephys/venus_constant')
#CST=np.loadtxt('filephys/Proxb_constant')
#CST=np.loadtxt('filephys/GJ504b_constant')
g=CST[0]
TT0=CST[3]
pref=CST[4]
cp=CST[1]
MM=CST[2]
R=8.314/(MM*1e-3)
gam=R/cp
nu=0.35
T0=460.
TMK=((Tm+TT0)**(nu)-(nu*T0**(nu)*np.log((pref/p)**gam)))**(1/nu)
TZK=((Tz+TT0)**(nu)-(nu*T0**(nu)*np.log((pref/p)**gam)))**(1/nu)
#deltaT=TZK-TMK

deltaW=[]
deltaT=[]
for i in range(len(Wxy)):
	deltaW.append(Wxy[i]-Wm)

for i in range(len(p)) :
        #print TZK[i],TMK[i]
        deltaT.append(TZK[i]-TMK[i])
#exit()
XX=np.linspace(0,60.4,len(Wxy))
pl.x = XX
pl.y = XX
pl.f = deltaW
pl.title = 'Vertical wind perturbation (m/s)'
pl.xlabel = 'x axis (km)'
pl.ylabel = 'y axis (km)'
#pl.colorbar='seismic'
pl.colorbar='RdBu_r'
pl.fmt="%.2f"
pl.vmin=-0.1
pl.vmax=0.1
pl.makeshow()

ppplot.figuref(x=10,y=8)
p1.x = XX
p1.f = deltaW[xx]
p1.ylabel = 'Vertical wind perturbation (K)'
p1.xlabel = 'y axis (km)'
p1.fmt="%.3f"
#pl.xmin=0
#pl.xmax=36
#p1.title = "Temperature perturbation at 3e4 Pa and x=2 km"
p1.makeshow()

ppplot.figuref(x=10,y=8)
p1.x = deltaT[x1:x2]
p1.f = p[x1:x2]
p1.ylabel = 'Vertical wind perturbation (K)'
p1.xlabel = 'y axis (km)'
p1.fmt="%.3f"
#pl.xmin=0
#pl.xmax=36
p1.logy = True
p1.invert= True
#p1.title = "Temperature perturbation at 3e4 Pa and x=2 km"
p1.makeshow()
#exit()
np.savetxt('deltat.txt',deltaT)
np.savetxt('deltaw.txt',deltaW[xx])
exit()
for i in range(11) : 
     np.savetxt('deltaw'+str(i),deltaW[xx+10*i])
#ifilename="./gif/W_"+str(i)
exit()
lamb=5 #(km)
N2=1e-4
wnl=(2*np.pi)/(lamb*1e3)
wnh2=2*(wnl**2)
wnv=(2*np.pi)/(1*1e3)
wn2=wnh2+wnv**2
omega2=(N2*wnh2)/wn2
calpha2=omega2/N2
calpha=sqrt(calpha2)
alpha=acos(calpha)*(180./np.pi)
print alpha
cgh=(sqrt(N2)*wnv)/(wn2)*(wnv)/sqrt(wn2)
cgz=(sqrt(N2)*wnv)/(wn2)*(-(sqrt(wnh2)))/sqrt(wn2)
print sqrt(omega2),alpha,cgh,cgz
