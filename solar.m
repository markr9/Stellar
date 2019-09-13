format compact
format long
close all
clear all

%solar model
Rsun=695700000 %[m]
Msun=1.989e30 %[kg]
rhomeansun=1408 %[kg/m^3]
Teffsun=5772 %[K]
Lsun=3.846e26 %[W]
X=0.73
Y=0.25
Z=0.02
kappa=0.02*(1+X) %[m^2kg^-1]
muion=1/(X+Y/4+Z/20)
mue=1/(1/2*(1+X))
mu=1/(1/muion+1/mue)

%constants
G=6.674e-11
c=299792458
u=1.66053904e-27
mh=1.67e-27
k=1.3806485e-23
o=5.670373e-8
a=(4*o)/c

%Lane-Emdun equation for n=3
% function lamed
z0=0, o0=1, o_dot=0 %o=x, z=y
zmax=7
n=3 %poltropic model
t=1e-9
% options=odeset('Events',@events);
f=@(z,o)[o(2);-(2/z*o(2)+o(1)^n)];
[z1,o1]=ode45(f,[t,zmax],[o0;z0]);
figure
plot(z1,o1(:,1))
grid on
xlabel('z')
ylabel('o')
title('Lane-Emden equation n=3')
% function dydt=f(z1,o1)
% dydt=[o(2);-(2/z*o(2)+o(1)^n)]
% function [value,isterminal,direction]=evenets(z1,o1)
% value=o1(1);
% isterminal=1;
% direction=-1;
%['NonNegative']

%polytropic
n=3
Rn=6.90
Dn=54.18
Mn=2.02
Bn=0.157

%solar test
rhoc=rhomeansun*Dn
alpha=Rsun/Rn
M=(4*pi*(Rsun/Rn)^3*Mn*rhomeansun*Dn)
M=(4*pi*alpha^3*Mn*rhomeansun*Dn)
Pc=(4*pi)^(1/3)*Bn*G*Msun^(2/3)*rhoc^(4/3)
Tc=(Pc*mu*u)/(rhoc*k)
Pgc=(rhoc*k*Tc)/(mu*u)
Prad=a/3*Tc^4
Pg=(rhoc*k*Tc)/(mu*u)
Pdeg=1.0036e7*(rhoc/mue)^(5/3)
Pdegrel=1.2435e10*(rhoc/mue)^(4/3)
PrPc=Prad/Pc
PrPg=Prad/(Pc+Prad)
PrPg2=(a*mu*u*Tc^3)/(3*rhoc*k)
Tc2=(2*G*u*Msun)/(3*Rsun*k) %3/2*k*T=(G*M*m)/R
Pc2=(G*Msun^2)/(4*pi*Rsun^4) %P=F/A: F=GM^2/R^1, A=4*pi*R^2
rho2=(Pc2*mu*u)/(k*Tc2)
beta2=1-PrPc
beta=1-PrPg

%lifetime
mfrac=0.15 %mass 1/2 L produced
eff=0.007
tau=((mfrac*eff*X*Msun*c^2)/Lsun)/(60^2*24*365.25)

%Lane-Emdun equation for n=3
figure
plot(z1,o1(:,1))
grid on
xlabel('z')
ylabel('o')
title('Lane-Emden equation n=3')
figure
plot(z1/Rn,o1(:,1))
grid on
xlabel('z')
ylabel('o')
title('Dimentionless n=3')

%profiles
rho1=rhoc*o1(:,1).^n;
P1=Pc*o1(:,1).^(n+1);
T1=Tc*o1(:,1);

figure
plot(z1/Rn,rho1/rhoc)
grid on
xlabel('z/Rn')
ylabel('rho/rhoc')
title('Denisty n=3')
figure
plot(z1*alpha/1000,rho1)
grid on
xlabel('r [km]')
ylabel('rho [kg/m^3]')
title('Denisty n=3')

figure
plot(z1/Rn,P1/Pc)
grid on
xlabel('z/Rn')
ylabel('P/Pc')
title('Pressure n=3')
figure
plot(z1*alpha/1000,P1)
grid on
xlabel('r [km]')
ylabel('P [Pa]')
title('Pressure n=3')

figure
plot(z1/Rn,T1/Tc)
grid on
xlabel('z/Rn')
ylabel('T/Tc')
title('Temperature n=3')
figure
plot(z1*alpha/1000,T1)
grid on
xlabel('r [km]')
ylabel('T [K]')
title('Temperature n=3')

%mass
r1=alpha*z1;
dr1=zeros(length(z1)-1,1); %84
i=1;
while i<length(z1)
    dr1(i)=r1(i+1)-r1(i);
    i=i+1;
end
j=1;
m1=zeros(length(z1)-2,1);
m2=0;
r2=zeros(length(z1)-2,1);
while j<(length(z1)-1)
    m11=4*pi*(rho1(j)+rho1(j+1))/2*((r1(j)+r1(j+1))/2)^2*dr1(j);
    r2(j)=(r1(j)+r1(j+1))/2;
    m2=m2+m11;
    m1(j)=m2;
    j=j+1;
end
m2
figure
plot(r2/alpha/Rn,m1/M)
grid on
xlabel('z/Rn')
ylabel('m/M')
title('Mass n=3')
figure
plot(r2/1000,m1)
grid on
xlabel('r [km]')
ylabel('m [kg]')
title('Mass n=3')

%mass2
dz1=zeros(length(z1)-1,1);
i=1;
while i<length(z1)
    dz1(i)=z1(i+1)-z1(i);
    i=i+1;
end
j=1;
m1=zeros(length(z1)-2,1);
m2=0;
z2=zeros(length(z1)-2,1);
while j<(length(z1)-1)
    m11=4*pi*alpha^3*rhoc*((z1(j)+z1(j+1))/2)^2*((o1(j)+o1(j+1))/2)^n*dz1(j);
    z2(j)=(z1(j)+z1(j+1))/2;
    m2=m2+m11;
    m1(j)=m2;
    j=j+1;
end
m2
figure
plot(z2/Rn,m1/M)
grid on
xlabel('z/Rn')
ylabel('m/M')
title('Mass n=3')
figure
plot(z2*alpha/1000,m1)
grid on
xlabel('r [km]')
ylabel('m [kg]')
title('Mass n=3')
m21=m2;

%mass3
dz1=zeros(length(z1)-1,1);
i=1;
while i<length(z1)
    dz1(i)=z1(i+1)-z1(i);
    i=i+1;
end
j=1;
m1=zeros(length(z1)-2,1);
m2=0;
z2=zeros(length(z1)-2,1);
while j<(length(z1)-1)
    m11=-4*pi*alpha^3*rhoc*((o1(j+1,2)*z1(j+1)^2)-(z1(j)^2*(o1(j,2))));
    z2(j)=(z1(j)+z1(j+1))/2;
    m2=m2+m11;
    m1(j)=m2;
    j=j+1;
end
m2
figure
plot(z2/Rn,m1/M)
grid on
xlabel('z/Rn')
ylabel('m/M')
title('Mass n=3')
figure
plot(z2*alpha/1000,m1)
grid on
xlabel('r [km]')
ylabel('m [kg]')
title('Mass n=3')
m22=m2;
m2=(m21+m22)/2

% m1=-4*pi*alpha^3*rhoc*[cumtrapz(z1.^2,o1(:,2).^n)]
% figure
% plot(z1/Rn,m1/m1(85))
% grid on
% xlabel('z/Rn')
% ylabel('m (normalised)')
% title('Mass n=3')
% sum(m1)

%energy
c1=1
q0=0.0012
lsun=4*pi*alpha^3*c1*q0*(X/0.5)^2*rhoc/rhoc*(Tc/Tc)^4*1e6
dz1=zeros(length(z1)-1,1);
i=1;
while i<length(z1)
    dz1(i)=z1(i+1)-z1(i);
    i=i+1;
end
j=1;
l1=zeros(length(z1)-2,1);
l2=0;
z2=zeros(length(z1)-2,1);
while j<(length(z1)-1)
    l11=lsun*((z1(j)+z1(j+1))/2)^2*((o1(j)+o1(j+1))/2)^(n+4)*dz1(j);
    z2(j)=(z1(j)+z1(j+1))/2;
    l2=l2+l11;
    l1(j)=l2;
    j=j+1;
end
l2

c2=1
q0a=9.7e-6
lsuna=4*pi*alpha^3*c2*q0a*(X/0.5)*(Z/0.02)*rhoc/rhoc*(Tc/Tc)^16*1e6
dz1=zeros(length(z1)-1,1);
i=1;
while i<length(z1)
    dz1(i)=z1(i+1)-z1(i);
    i=i+1;
end
j=1;
l1a=zeros(length(z1)-2,1);
l2a=0;
z2=zeros(length(z1)-2,1);
while j<(length(z1)-1)
    l11a=lsuna*((z1(j)+z1(j+1))/2)^2*((o1(j)+o1(j+1))/2)^(n+16)*dz1(j);
    z2(j)=(z1(j)+z1(j+1))/2;
    l2a=l2a+l11a;
    l1a(j)=l2a;
    j=j+1;
end
l2a

figure
plot(z2/Rn,(l1+l1a)/Lsun)
grid on
xlabel('z/Rn')
ylabel('l/L')
title('Luminosity n=3')
figure
plot(z2*alpha/1000,l1+l1a)
grid on
xlabel('r [km]')
ylabel('l [W]')
title('Luminosity n=3')

%energy2
l=zeros(length(z1),1);
i=1;
while i<length(z1)+1
    l(i)=-(4*a*c*Tc^4*4*pi*alpha)/(3*kappa*rhoc)*z1(i)^2*o1(i,2);
    i=i+1;
end
figure
plot(z1/Rn,l/Lsun)
grid on
xlabel('z/Rn')
ylabel('l/L')
title('Luminosity n=3')
figure
plot(z1*alpha/1000,l)
grid on
xlabel('r [km]')
ylabel('l [W]')
title('Luminosity n=3')

%eddington model
kk=0.003
beta1=0.25:1e-6:1;

x=sqrt((1-beta1)./(kk*beta1.^4));
%x=mu^2*M/Msun

figure
plot(log10(x),beta1)
xlabel('log(mu^2*M/Msun)')
ylabel('Beta')
grid on

figure
Medd=x/mu^2;
plot(log10(Medd),beta1)
xlabel('log(M/Msun)')
ylabel('Beta')
grid on

betasun=0.9996041
Ledd=(4*pi*c*G*Msun)/kappa
L=Ledd*(1-betasun)
L=L/Lsun

Pcedd=((3*(1-betasun))/a)^(1/3)*(k/(betasun*mu*u))^(4/3)*rhoc^(4/3)
Tcedd=(betasun*Pcedd*mu*u)/(rhoc*k)

%K
K=Pc/rhoc^(4/3)
K2=(M/(4*pi*Mn))^(2/3)*pi*G
K3=((3*(1-betasun))/a)^(1/3)*(k/(betasun*mu*u))^(4/3)
K4=(alpha^2*4*pi*G*rhoc^(2/3))/(n+1)

%Teff
Teff=(Lsun/(4*pi*alpha^2*Rn^2*o))^(1/4)
Teff2=(Lsun/(4*pi*Rsun^2*o))^(1/4)