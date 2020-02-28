function [pass,maxerr] = test(opt)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 80;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.3]);
rng(2)
K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.2);

alpha = 0.0005;

if opt.Display
    figure(8),clf
    axhandle = axes();
else
    axhandle = [];
end

Pfit1 = obir(S,K,r,'tikh',alpha,'NoiseLevelAim',0.05,'TolFun',1e2,'MaxFunEvals',100,'MaxOuterIter',100,'MaxIter',100,'RegOrder',1,'HuberParam',15);
Pfit2 = obir(S,K,r,'tikh',alpha,'NoiseLevelAim',0.05,'TolFun',1e1,'MaxFunEvals',500,'MaxOuterIter',500,'MaxIter',5000,'RegOrder',2,'HuberParam',1.35);

err = norm(P - Pfit1) < norm(P - Pfit2);
maxerr = max(abs(Pfit1 - Pfit2));
 

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Pfit,'b')
    legend('truth','OBIR','Tikh')
end

end