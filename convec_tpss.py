#! /usr/bin/env python


import numpy as np
from ppclass import pp
import ppplot
import ppcompute
from math import *
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy import signal
#from planets import Venus

vv="W"
tt=30
xx=90
xm="0,1e10"
yy=90
zz=80
tp="T"
#fi="Venus4070_CKJan15H_0_minuit_0/wrfout_d01_9999-01-01_00:00:00"
#fi="Venus4070_CKJan15H_0_minuit_0/wrfout_d01_9999-01-01_00:00:00"
#fi="Venus4070_CKJan15H_0_minuit_0/wrfout_d01_9999-01-01_22:13:20"
#fi="Venus4070_CKMai14Z_0_minuit/wrfout_d01_9999-01-02_03:46:40"
#fi="Venus4070_1d_0_minuit/wrfout_d01_9999-01-01_00:00:00"
#fi="Venus4070_CKJan15H_0_midi_100tes/wrfout_d01_9999-01-01_22:13:20"
#fi="Venus090_phys_0_midi5/wrfout_d01_9999-01-01_00:00:00"
#fi='Venusphtest_hr_cpcst2/wrfout_d01_9999-01-01_00:00:00'
#fi='Venusphtest_hr_cpvar2/wrfout_d01_9999-01-01_00:00:00'
#fi='GJ504bTest5/wrfout_d01_9999-01-01_05:33:20'
#fi='Test_dyncpscte_physcpvar2/wrfout_d01_9999-01-01_00:00:00'
#fi="V_eq_n_ws/wrfout_d01_9999-01-01_00:00:00"
#fi="GJ504b/wrfout_d01_9999-01-01_11:06:40"
fi="GJ504b2/wrfout_d01_9999-01-01_16:40:00"
#fi="wrfout_d01_9999-01-03_07:33:20"
#fi="VENUS_PBL5/wrfout_d01_9999-01-01_00:00:00"
fi="wrfout_d01_9999-01-01_22:13:20"

p1 = ppplot.plot1d()
p11 = ppplot.plot1d()
pl = ppplot.plot2d()
pl2 = ppplot.plot2d()

T,x,y,z,t = pp(file=fi,var=tp,t=tt,x=xm,y=xm).getfd()
#T2 = pp(file=fi,var=tp,t=95,x=xm,y=xm).getf()

p = pp(file=fi,var="PTOT",t=tt,x=xm,y=xm).getf()
#lw = pp(file=fi,var="HR_LW",t=tt,x=xm,y=xm).getf() 
#dtrad = pp(file=fj,var="HR_LW",t=tt,x=xm,y=xm).getf()
#hr = open("input_hr", "w")
CST=np.loadtxt('filephys/venus_constant')
#CST=np.loadtxt('filephys/Proxb_constant')
#CST=np.loadtxt('filephys/GJ504b_constant')
g=CST[0]
TT0=CST[3]
pref=CST[4]
cp=CST[1]
MM=CST[2]
R=8.314/(MM*1e-3)
#gam=R/cp
#TT0=735
#TT0=1000.
#pref = 92.e5
#pref = 6.8e6
gam=192./1000.
#gam1=192./900.
teta=T+TT0
#teta2=T2+600.
nu=0.35
T0=460.
#g=8.87
#g=100.
#R=191.84383904727036
#R=8.314/(2.35*1e-3)
#cp=1.1787E4
#cp=1000.*((temp/T0)**nu)
gam1=R/cp
print R,gam1

temp=(teta**(nu)-(nu*T0**(nu)*np.log((pref/p)**gam)))**(1/nu)
#temp2=(teta2**(nu)-(nu*T0**(nu)*np.log((pref/p)**gam)))**(1/nu)
temp1=teta*(p/pref)**(gam1)
rho=p/(191.84*temp1)
#cp=1000.*((temp/T0)**nu)
#cp=1.1787E4
#print cp
Tteta=pp(file=fi,var="T",y=yy,t=tt).getf()
TTT=pp(file=fi,var="T",t=tt).getf()
MM=[]
TT=[]
TTTK=[]
pertTTTK=[]
DTM=[]
DT=[]
DT2=[]
TTT0 = pp(file=fi,var=tp,t=93,x=xx,y=xx).getf()
teta0=TTT0+TT0
temp0=(teta0**(nu)-(nu*T0**(nu)*np.log((pref/p)**gam)))**(1/nu)
for i in range(6) :
    Ti = pp(file=fi,var=tp,t=93+i,x=xx,y=xx).getf()
    Ti2 = pp(file=fi,var=tp,t=93+i+1,x=xx,y=xx).getf()
    tetai=Ti+TT0
    tetai2=Ti2+TT0
    tempi=(tetai**(nu)-(nu*T0**(nu)*np.log((pref/p)**gam)))**(1/nu)
    tempi2=(tetai2**(nu)-(nu*T0**(nu)*np.log((pref/p)**gam)))**(1/nu)
    DTM.append(np.amax(tempi2-tempi)-np.amin(tempi2-tempi))
    DT.append(tempi2-tempi)
    DT2.append(tempi-temp0)
    #DT.append(np.amax(tempi)-np.amin(tempi))
    #print ppcompute.max(tempi2-tempi)
    #print ppcompute.min(tempi2-tempi)
#exit()
TEMPS=[3.3,6.6,10,13.3,15.6,19]

p1.f = DTM
p1.x = TEMPS
p1.xlabel = 'Time (min)'
p1.ylabel = 'Tmax-Tmin (K)'
p1.fmt="%.2f"
#p1.logy = True
#p1.title = 'Temperature variation over 15 min (40-70 km) No convection'
#p1.invert = True
p1.makeshow()

pl2.f = DT2
pl2.y = p
pl2.x = TEMPS
pl2.ylabel = 'Pressure (Pa)'
pl2.xlabel = 'Time (min)'
pl2.logy = True
pl2.invert = True
#pl2.title = "Vertical wind pertubation (m/s)"
pl2.vmin = -4E-1
pl2.vmax = 4E-1
pl2.colorbar='RdBu_r'
pl2.makeshow()
exit()
for i in range(len(p)):
	#TT.append((np.array(Tteta[i]+TT0)**(nu)-(nu*T0**(nu)*np.log((pref/p[i])**gam)))**(1/nu))
        TT.append((Tteta[i]+TT0)*(p[i]/pref)**(gam1))
        #for j in range(len(Tteta[1])) :
        TTTK.append((np.array(TTT[i]+TT0)**(nu)-(nu*T0**(nu)*np.log((pref/p[i])**gam)))**(1/nu))
        pertTTTK.append(TTTK[i]-temp[i])

alpha=0.025/(rho*cp)
PP = pp(file=fi,var="PTOT",t=tt,y=yy).getf()
Pm_e = pp(file=fi,var="PTOT",t=tt,x=xm,y=xm).getf()
Tm_e = pp(file=fi,var=tp,t=tt,x=xm,y=xm).getf()

#print p, temp, teta
#exit()
#file1=open("temp.txt","w")
#file2=open("temp2.txt","w")
#file1.write(temp1)
#file2.write(temp2)
#file1.close()
#file2.close()
#print p
p1.f = p
p1.x = temp1
p1.xlabel = 'Spatial mean temperature (K)'
p1.ylabel = 'Pressure (Pa)'
p1.logy = True
#p1.title = 'Temperature variation over 15 min (40-70 km) No convection'
p1.invert = True
p1.makeshow()
p1.f = p
p1.x = teta 
p1.xlabel = 'Spatial mean potential temperature (K)'
p1.ylabel = 'Pressure (Pa)'
p1.logy = True
#p1.title = 'Potential Temperature'
#p1.xmin=720
#p1.xmax=1080
p1.invert = True
p1.makeshow()
#dt=ppcompute.deriv1d(tetap,z)

dtmdp=ppcompute.deriv1d(temp,p)
dtdz=(-p*g/(191.84383904727036*temp))*dtmdp
dtmdp1=ppcompute.deriv1d(temp1,p)
dtdz1=(-p*g/(R*temp1))*dtmdp1
#Ss=(temp/tetap)*dt
ssm=dtdz+(g/cp)
ssm1=dtdz1+(g/cp)
#Ss=dt+(8.87/900.)

p1.f = p
p1.x = 1000*ssm1
p1.xlabel = 'Static Stability (K/km)'
p1.ylabel = 'Pressure (Pa)'
p1.logy = True
p1.invert = True
p1.makeshow()

#p1.f = p
#p1.x = dtrad*100
#p1.xlabel = 'dTrad (K/s)'
#p1.ylabel = 'Pressure (Pa)'
#p1.logy = True
#p1.invert = True
#p1.makeshow()
#exit()
#dtetadp=ppcompute.deriv1d(Tm_e,p)
#dtetadz=(-p*8.87/(191.84383904727036*temp))*dtetadp
#N2=(8.87/teta)*dtetadz
N2=(g/temp)*ssm
N21=(g/temp1)*ssm1
#print len(N2),len(p), len(z)
#exit()
p1.f = p 
p1.x = 10000*N2
p1.ylabel = 'Pressure (Pa)'
p1.xlabel = 'Square of the mean Brunt-Vaisala frequency (1e-4 s-2)'
p1.logy = True
p1.invert = True
p1.makeshow()

pertT=[]
pertP=[]
vmf=[] 	
mmf=[]
chfu=[]
chfw=[]
vef=[]
baker=[]
vmfi=[]
mmfi=[]
chfui=[]
chfwi=[]
vefi=[]
EPf=[]

#ro=p/(temp1*191.84383904727036)
ro=p/(temp1*R)
WWW=pp(file=fi,var="W",t=tt).getf()
w,x,y,zw,t = pp(file=fi,var="W",t=tt,y=yy).getfd()
wz = pp(file=fi,var="W",t=tt,z=250).getf()
wp_e = pp(file=fi,var="W",t=tt,x=xm,y=xm).getf()
u,x,y,z,t = pp(file=fi,var="U",t=tt,y=yy).getfd()
up_e = pp(file=fi,var="U",t=tt,x=xm,y=xm).getf()
v,x,y,z,t = pp(file=fi,var="V",t=tt,y=yy).getfd()
vp_e = pp(file=fi,var="V",t=tt,x=xm,y=xm).getf()
pertu=[]
pertv=[]
pertw=[]
pertwz=[]
perttp=[]
MAXtp=[]
MINtp=[]
MAXw=[]
MINw=[]
#print len(w),len(u),len(v),len(wp_e[0]),len(up_e[0]),len(vp_e[0])
tteta=pp(file=fi,var=tp,t=tt,x=xx,y=yy).getf()
#TK=((tteta+600)**(nu)-(nu*T0**(nu)*np.log((pref/p)**gam)))**(1/nu)
TK=(tteta+TT0)*((p/pref)**gam1)
for i in range(len(p)):
	pertu.append(u[i]-up_e[i])
	pertv.append(v[i]-vp_e[i])
	pertw.append(w[i]-wp_e[i])
        #pertwz.append(wz[i]-wp_e[50])
        pertT.append(TT[i]-temp1[i])
        perttp.append(TK[i]-temp1[i])
        MAXtp.append(np.amax(pertTTTK[i]))
        MAXw.append(np.amax(WWW[i]))
	MINtp.append(np.amin(pertTTTK[i]))
        MINw.append(np.amin(WWW[i]))
	pertP.append(PP[i]-Pm_e[i])
        EPf.append(-1*np.mean(np.array(pertw[i])*(np.array(pertu))[i,1:]))
        #vmf.append(np.mean(rho[i]*np.array(pertw[i])*(np.array(pertu))[i,1:]))
        #mmf.append(np.mean(pertw[i]*pertv[i]))
        #chfu.append(np.mean(pertT[i]*(np.array(pertu))[i,1:]))
	chfw.append(rho[i]*cp*np.mean(pertT[i]*np.array(pertw)[i]))
	#vef.append(np.mean(pertP[i]*np.array(pertw)[i]))
	#baker.append(rho[i]**0.5*np.array(pertw)[i])
	#vmfi.append(rho[i]*np.array(pertw)[i]*(np.array(pertu))[i,1:])
	#mmfi.append(np.array(pertw)[i]*pertv[i])
	#chfui.append(pertT[i]*(np.array(pertu))[i,1:])
	#chfwi.append(pertT[i]*np.array(pertw)[i])
	#vefi.append(pertP[i]*np.array(pertw)[i])
print len(np.array(pertu)),len(np.array(pertu)[0])
dz=ppcompute.deriv1d(chfw,p)
dzchfw=(-p*8.87/(191.84383904727036*temp1))*dz
X=np.linspace(0,1800,len(x))

p1.f = p
p1.x = np.array(MAXtp)
p1.ylabel = 'Pressure (Pa)'
p1.xlabel = 'K'
p1.title = 'Temperature pertubation profil'
p1.logy = True
p1.invert = True
p1.makeshow()		

p1.f = p
p1.x = np.array(MINtp)
p1.ylabel = 'Pressure (Pa)'
p1.xlabel = 'K'
p1.title = 'Temperature pertubation profil'
p1.logy = True
p1.invert = True
p1.makeshow()

#p1.f = p
#p1.x = np.array(perttp)
#p1.ylabel = 'Pressure (Pa)'
#p1.xlabel = 'Temperature pertubation profil (K)'
#p1.title = 'Temperature pertubation profil'
#p1.logy = True
#p1.invert = True
#p1.makeshow()
#np.savetxt('deltaTvert',perttp)

#pl.f = np.array(pertT)
#pl.y = z
#pl.x = X
#pl.ylabel = 'Pressure (Pa)'
#pl.xlabel = 'x (km)'
#pl.title = "Temperature pertubation (K)"
#pl.logy = True
#pl.invert = True
#pl.vmin = -1E-1
#pl.vmax = 1E-1
#pl.colorbar='RdBu_r'
#pl.makeshow()

#exit()
#pl.f = np.array(w[1:])
#pl.y = z
#pl.x = X
#pl.ylabel = 'Pressure (Pa)'
#pl.xlabel = 'x (km)'
#pl.title = "Vertical wind (m/s)"
#pl.logy = True
#p#l.invert = True
#p#l.vmin = -1E-1
#pl.vmax = 1E-1
#pl.colorbar='RdBu_r'
#pl.makeshow()

#pl2.f = np.array(pertwz)
#pl2.y = X
#pl2.x = X
#pl2.ylabel = 'y (km)'
#pl2.xlabel = 'x (km)'
#pl2.title = "Vertical wind pertubation (m/s)"
#pl.vmin = -1E-1
#pl.vmax = 1E-1
#pl2.colorbar='RdBu_r'
#pl2.makeshow()
#print EPf
#print z
p1.x = EPf
p1.y = z
p1.ylabel = 'Pressure (Pa)'
p1.xlabel = "EP flux"
#p1.title = "VMF rho.w'.u'"
p1.logy = True
p1.invert = True
#p1.vmin = -5E-1
#p1.vmax = 5E-1
p1.makeshow()
#exit()
p1.x = vmf
p1.y = z
p1.ylabel = 'altitude p'
p1.xlabel = "VMF rho.w'.u'"
#p1.title = "VMF rho.w'.u'"
p1.logy = True
p1.invert = True
#p1.vmin = -5E-1
#p1.vmax = 5E-1
p1.makeshow()

p1.x = mmf
p1.y = z
p1.ylabel = 'altitude p'
p1.xlabel = "MMF w'v'"
#p1.title = "MMF w'v'"
p1.logy = True
p1.invert = True
#p1.vmin = -5E-1
#p1.vmax = 5E-1
p1.makeshow()

p1.x = chfu
p1.y = z
p1.ylabel = 'altitude p'
p1.xlabel = "chfu T'u'"
#p1.title = "chfu T'u'" 
p1.logy = True
p1.invert = True
#p1.vmin = -5E-1
#p1.vmax = 5E-1
p1.makeshow()

p1.x = chfw
p1.y = z
p1.ylabel = 'Pressure (Pa)'
p1.xlabel = 'Convective heat flux (W/m2)'
#p1.title = "Turbulent heat flux"
p1.logy = True
p1.invert = True
#p1.vmin = -5E-1
#p1.vmax = 5E-1
p1.makeshow()

p1.x = dzchfw
p1.y = z
p1.ylabel = 'Pressure (Pa)'
p1.xlabel = 'Grandient of Convective heat flux (W/m2/m)'
#p1.title = "Turbulent heat flux"
p1.logy = True
p1.invert = True
#p1.vmin = -5E-1
#p1.vmax = 5E-1
p1.makeshow()

p1.x = vef
p1.y = z
p1.ylabel = 'altitude p'
p1.xlabel = 'x (km)'
p1.title = "chfw p'w'"
p1.logy = True
p1.invert = True
#p1.vmin = -5E-1
#p1.vmax = 5E-1
p1.makeshow()

p1.x = vef
p1.y = z
p1.ylabel = 'altitude p'
p1.xlabel = 'baker'
p1.title = "rho^0.5 * w'"
p1.logy = True
p1.invert = True
#p1.vmin = -5E-1
#p1.vmax = 5E-1
p1.makeshow()
exit()
pl.f = vmfi
pl.x = np.linspace(0,36,180)
pl.y = z
pl.ylabel = 'altitude p'
pl.xlabel = 'x (km)'
pl.title = "w'u'"
pl.logy = True
pl.invert = True
#pl.vmin = -5E-1
#pl.vmax = 5E-1
pl.makeshow()

pl.f = mmfi
pl.x = np.linspace(0,36,180)
pl.y = z
pl.ylabel = 'altitude p'
pl.xlabel = 'x (km)'
pl.title = "w'v'"
pl.logy = True
pl.invert = True
#pl.vmin = -5E-1
#pl.vmax = 5E-1
pl.makeshow()

pl.f = chfui
pl.x = np.linspace(0,36,180)
pl.y = z
pl.ylabel = 'altitude p'
pl.xlabel = 'x (km)'
pl.title = "chfu T'u'" 
pl.logy = True
pl.invert = True
#pl.vmin = -5E-1
#pl.vmax = 5E-1
pl.makeshow()

pl.f = chfwi
pl.x = np.linspace(0,36,180)
pl.y = z
pl.ylabel = 'altitude p'
pl.xlabel = 'x (km)'
pl.title = "chfw T'p' "
pl.logy = True
pl.invert = True
#pl.vmin = -5E-1
#pl.vmax = 5E-1
pl.makeshow()

pl.f = vefi
pl.x = np.linspace(0,36,180)
pl.y = z
pl.ylabel = 'altitude p'
pl.xlabel = 'x (km)'
pl.title = "chfw p'w'"
pl.logy = True
pl.invert = True
#pl.vmin = -5E-1
#pl.vmax = 5E-1
pl.makeshow()


