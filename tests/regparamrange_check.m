function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);

K = dipolarkernel(t,r);
L = regoperator(Dimension,2);
logRes = 0.1;
alpha = regparamrange(K,L,'logResolution',logRes);

err = (length(alpha)~=85);
maxerr = NaN;
data = [];



end