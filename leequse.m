%Lane-Emdun equation for n=3
z0=0, o0=1, o_dot=0 %o=x, z=y
zmax=8
n=3
t=1e-99
f=@(z,o)[o(2);-(2/z*o(2)+o(1)^n)];
[z1,o1]=ode45(f,[t,zmax],[o0;z0],[]);
plot(z1,o1(:,1))
grid on
xlabel('z/R')
ylabel('o/pc')
title('Lane-Emden equation n=3')