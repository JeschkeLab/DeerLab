function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 80;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
P = P/sum(P);

K = dipolarkernel(t,r);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = K*P;

[OptParam,Functionals,RegParams] = selregparam(DipEvoFcn,K,RegMatrix,'tikhonov',{'aic','gcv','lr'});
%Accept testif all values are the same (should be as there is no noise)
err(1) = any(diff(OptParam) > 1e-2);
err(2) = any(abs(OptParam - 0.002) > 1e-4);
err = any(err);
data = [];
maxerr = max(abs(OptParam - 0.002));

if opt.Display
   figure(8),clf
   hold on
   plot(RegParams,Functionals{1}/max(Functionals{1}),'.')
   plot(RegParams,Functionals{2}/max(Functionals{2}),'.')
end


end