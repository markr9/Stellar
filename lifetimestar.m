format compact

%solar model
Rsun=695700000 %[m]
Msun=2.1554e30 %[kg]
rhomeansun=1408 %[kg/m^3]
Teffsun=5975 %[K]
Lsun=5.92e26 %[W]
X=0.73
Y=0.25
Z=0.02
c=299792458

%lifetime
mfrac=0.15 %mass 1/2 L produced
format long
mH=1.00782503207
mHe=4.00260325415
eff=(4*mH-mHe)/(4*mH)
eff2=0.007
format short
tau=((mfrac*eff*X*Msun*c^2)/Lsun)/(60^2*24*365.25)
tau/2