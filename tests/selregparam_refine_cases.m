function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

N = 60;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.01);

alphaopt1 = selregparam(S,K,r,'tikhonov','aic','Range',linspace(0.1,0.4,20));
alphaopt2 = selregparam(S,K,r,'tikhonov','aic','Range',linspace(0.4,0.6,20),'Refine',true);
alphaopt3 = selregparam(S,K,r,'tikhonov','aic','Range',linspace(0.01,0.1,20),'Refine',true);


err(1) = any(abs(alphaopt1 - alphaopt2) > 3e-1);
err(2) = any(abs(alphaopt1 - alphaopt3) > 3e-1);
err(3) = any(abs(alphaopt2 - alphaopt3) > 3e-1);

err = any(err);
data = [];
maxerr = 0;

end