function [pass,maxerr] = test(opt)

% Check that obir() runs with specified options

rng(1)
t = linspace(0,3,200);
r = linspace(0,5,100);
P = dd_twogauss(r,[2,0.3,3.5,0.3,0.5]);
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.05);
alpha = 0.05;

Pfit1 = obir(S,K,r,'tikh',alpha,'NoiseLevelAim',0.05,'TolFun',1e-1,'MaxFunEvals',100,'MaxOuterIter',100,'MaxIter',100,'RegOrder',2,'HuberParam',1.4);
Pfit2 = obir(S,K,r,'tikh',alpha,'NoiseLevelAim',0.05,'TolFun',1e-2,'MaxFunEvals',500,'MaxOuterIter',500,'MaxIter',5000,'RegOrder',2,'HuberParam',1.35);

% Pass: the fit with the better options is better
pass = mean(abs(P - Pfit2)) < mean(abs(P - Pfit1));

maxerr = max(abs(P - Pfit2));
 
if opt.Display
    plot(r,P,'k',r,Pfit1,r,Pfit2)
    legend('truth','bad options','good options')
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    grid on, axis tight, box on
end

end