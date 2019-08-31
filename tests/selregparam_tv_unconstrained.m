function [err,data] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);
Distribution = Distribution/sum(Distribution);

K = dipolarkernel(t,r);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = K*Distribution;

RegParamSet = regparamrange(K,RegMatrix);
[OptParam,Functionals,RegParams] = selregparam(RegParamSet,DipEvoFcn,K,RegMatrix,'tv',{'aic','gcv'},'NonNegConstrained',false);

%Accept testif all values are the same (should be as there is no noise)
err = any(any(OptParam - OptParam' > 1e-2));
data = [];

if opt.Display
   figure(8),clf
   hold on
   plot(RegParamSet,Functionals{1})
   plot(RegParamSet,Functionals{2})
end


end