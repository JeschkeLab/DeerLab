function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.3]);
rng(2)
K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.2);

alpha = 0.05;

if opt.Display
    figure(8),clf
    axhandle = axes();
else
    axhandle = [];
end

Pfit1 = obir(S,K,r,'huber',alpha,'NoiseLevelAim',0.05,'RegOrder',1);
Pfit2 = obir(S,K,r,'huber',alpha,'NoiseLevelAim',0.05,'RegOrder',1,'HuberParameter',1.35);
err(1) = any(abs(Pfit1 - Pfit2)>1e-10);

Pfit1 = obir(S,K,r,'tikh',alpha,'NoiseLevelAim',0.05,'TolFun',1e2,'MaxFunEvals',100,'MaxOuterIter',100,'MaxIter',100);
Pfit2 = obir(S,K,r,'tikh',alpha,'NoiseLevelAim',0.05,'TolFun',1e-2,'MaxFunEvals',500000,'MaxOuterIter',5000,'MaxIter',500000);
err(2) = norm(P - Pfit1) < norm(P - Pfit2);


err = any(err);
maxerr = max(abs(Pfit1 - Pfit2));
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Pfit,'b')
    legend('truth','OBIR','Tikh')
end

end