function [pass,maxerr] = test(opt)

% Test time-domain global fit of parametric models

rng(2)
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
info = dd_twogauss();
mymodel = @(t,param)dipolarkernel(t,r)*dd_twogauss(r,param);
InitialGuess = [2 0.1 5 0.1 0.5];
range = [info.parameters(:).range];
upper = range(2:2:end);
lower = range(1:2:end-1);
parglobal = fitparamodel(Ss,mymodel,{t1,t2,t3},InitialGuess,'Upper',upper,'lower',lower);
parlocal = fitparamodel(S3,mymodel,t3,InitialGuess,'Upper',upper,'lower',lower);

Pglobal = dd_twogauss(r,parglobal);
Plocal = dd_twogauss(r,parlocal);


normResult = norm(P - Pglobal);
norm3 = norm(P - Plocal);


% Pass: global fit is a better fit than the local one
pass = all(normResult < norm3);
 
maxerr = max(P - Pglobal);

if opt.Display
    subplot(121)
    hold on
    plot(t1,S1,'k')
    plot(t1,K1*Pglobal,'r')
    plot(t2,S2+1,'k')
    plot(t2,K2*Pglobal + 1,'r')
    plot(t3,S3+2,'k')
    plot(t3,K3*Pglobal + 2,'r')
    grid on, axis tight, box on
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    subplot(122)
    hold on
    plot(r,P,'k')
    plot(r,Pglobal,'r')
    plot(r,Plocal,'b--')
    grid on, axis tight, box on
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
end

end