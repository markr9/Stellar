%HR dig.
R=1 %[Rsun]
T=52000:-1:2200; %[K]
L=R^2*(T/5778).^4;
loglog(T,L,'b-')
hold on
R=[0.0001,0.001,0.01,10,100,1000]
i=1;
while i<7
    L=R(i)^2*(T/5778).^4;
    loglog(T,L,'cy-')
    i=i+1;
end
xlabel('Temperature (log(K/Tsun)')
ylabel('Luminosity log(L/Lsun))')
title('HR digram')
grid on
hold off
hold on
plot(5778,1,'o')
hold off