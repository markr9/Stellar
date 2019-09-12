%Lane-Emdun equation for n=3
function [n,z1,o1,te,ye,ie]=lamed()
z0=0, o0=1, o_dot=0 %o=x, z=y
zmax=7
n=3 %polytropic model
t=1e-9
options=odeset('Events',@events);
[z1,o1,te,ye,ie]=ode45(@f,[t,zmax],[o0;z0],options);
te
ye
ie
figure
plot(z1,o1(:,1),te,ye(:,1),'r-o')
grid on
xlabel('z')
ylabel('o')
title('Lane-Emden equation n=3')
end
function dydt=f(z1,o1)
n=3;
dydt=[o1(2);-(2/z1*o1(2)+o1(1)^n)];
%f=@(z,o)[o(2);-(2/z*o(2)+o(1)^n)];
end
function [value,isterminal,direction]=events(z1,o1)
value=o1(1);
isterminal=1;
direction=-1;
end