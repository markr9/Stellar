%Pressure dominance
X=0.73
Y=0.25
Z=0.02
muion=1/(X+Y/4+Z/20)
mue=1/(1/2*(1+X))
mu=1/(1/muion+1/mue)
n=3
k=0.003
beta=0.25:1e-6:1;

x=sqrt((1-beta)./(k*beta.^4));
%x=mu^2*M/Msun

figure(1)
plot(log10(x),beta)
xlabel('log(mu^2*M/Msun)')
ylabel('Beta')
grid on

figure(2)
M=x/mu^2;
plot(log10(M),beta)
xlabel('log(M/Msun)')
ylabel('Beta')
grid on

