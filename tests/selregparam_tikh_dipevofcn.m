function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 80;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;

RegParamSet = regparamrange(Kernel,RegMatrix);
[OptParam,Functionals,RegParams] = selregparam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,'tikhonov',{'aic','gcv','lr'});
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