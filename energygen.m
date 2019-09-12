%pp and CNO
figure(1)
X=0.73
Xcno=0.01
rho=54.18*1408
T=0e6:0.001e6:10e7;
e0pp=2.6e-37
e0cno=7.9e-118
epp=e0pp*rho*X*X*T.^4.5;
ecno=e0cno*rho*X*Xcno*T.^16;
plot(T,log10(epp))
hold on
plot(T,log10(ecno))
hold off
xlabel('T (K)')
ylabel('Energy generation per kg per s')
grid on

%pp and CNO (better)
figure(2)
X=0.73
Xcno=0.01
rho=54.18*1408
T=0:0.001:100;
e0pp=1.08e-12
e0cno=8.24e-31
epp=e0pp*rho*X*X*T.^4;
ecno=e0cno*rho*X*Xcno*T.^19.9;
plot(T*10^6,log10(epp))
hold on
plot(T*10^6,log10(ecno))
hold off
xlabel('T (K)')
ylabel('Energy generation per kg per s')
grid on

%lecture intenet
figure(3)
X=0.73
Xcno=0.01
rho=54.18*1408
T=0e6:0.001e6:100e6;
epp=0.0012*(X/0.5)^2*(rho/10^5)*(T/1.5e7).^4;
ecno=9.7e-6*(X/0.5)*(Xcno/0.0029)*(rho/10^5)*(T/1.5e7).^16;
plot(T,log10(epp))
hold on
plot(T,log10(ecno))
hold off
xlabel('T (K)')
ylabel('Energy generation per kg per s')
grid on

figure(4)
ratio=epp./ecno;
plot(T,log(ratio))
xlabel('T (K)')
ylabel('log(epp/ecno)')
grid on
