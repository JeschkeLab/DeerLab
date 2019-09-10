function [err,data] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 80;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = K*P;

RegParamSet = regparamrange(K,RegMatrix);
[OptParam,~,~] = selregparam(DipEvoFcn,K,RegMatrix,'huber',{'aic','gcv'});

%Accept testif all values are the same (should be as there is no noise)
err = any(any(OptParam - OptParam' > 1e-2));
data = [];



end