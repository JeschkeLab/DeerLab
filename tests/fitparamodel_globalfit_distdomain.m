function [pass,maxerr] = test(opt)

% Test distance-domain global fit of parametric models

dt = 0.008;
r = linspace(1,5,100);
P = dd_twogauss(r,[2,0.3,4,0.3,0.5]);

Ntime1 = 100;
t1 = linspace(0,dt*Ntime1,Ntime1);
K1 = dipolarkernel(t1,r);
S1 = K1*P + whitegaussnoise(Ntime1,0.03);

Ntime2 = 200;
t2 = linspace(0,dt*Ntime2,Ntime2);
K2 = dipolarkernel(t2,r);
S2 = K2*P + whitegaussnoise(Ntime2,0.05);

Ntime3 = 300;
t3 = linspace(0,dt*Ntime3,Ntime3);
K3 = dipolarkernel(t3,r);
S3 = K3*P + whitegaussnoise(Ntime3,0.1);

Ss = {S1,S2,S3};
Ks = {K1,K2,K3};

[~,Pglobal] = fitparamodel(Ss,@dd_twogauss,r,Ks);
[~,Plocal1] = fitparamodel(S1,@dd_twogauss,r,K1);
[~,Plocal2] = fitparamodel(S2,@dd_twogauss,r,K2);
[~,Plocal3] = fitparamodel(S3,@dd_twogauss,r,K3);

rmsdglobal = norm(P - Pglobal);
rmsdlocal1 = norm(P - Plocal1);
rmsdlocal2 = norm(P - Plocal2);
rmsdlocal3 = norm(P - Plocal3);

% Pass: global fit is a better fit than the local one
pass = all(rmsdglobal < [rmsdlocal1 rmsdlocal2 rmsdlocal3]);
 
maxerr = max(abs(P - Pglobal));

if opt.Display
subplot(121)
hold on
plot(t1,S1,'k')
plot(t2,S2+1,'k')
plot(t3,S3+2,'k')
plot(t1,K1*Pglobal,'r')
plot(t2,K2*Pglobal + 1,'r')
plot(t3,K3*Pglobal + 2,'r')
grid on, axis tight, box on
xlabel('r [nm]')
ylabel('P(r) [nm^{-1}]')
subplot(122)
hold on
plot(r,P,'k')
plot(r,Pglobal,'r')
plot(r,Plocal1,'g--')
plot(r,Plocal2,'b--')
plot(r,Plocal3,'r--')
grid on, axis tight, box on
xlabel('r [nm]')
ylabel('P(r) [nm^{-1}]')
end

end