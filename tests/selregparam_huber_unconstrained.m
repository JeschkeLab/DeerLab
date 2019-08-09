function [err,data] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;

RegParamSet = regparamrange(Kernel,RegMatrix);
[OptParam,Functionals,RegParams] = selregparam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,{'aic','gcv'},'RegType','huber','NonNegConstrained',false);

%Accept testif all values are the same (should be as there is no noise)
err = any(any(OptParam - OptParam' > 1e-2));
data = [];

if opt.Display
   figure(8),clf
   hold on
   plot(RegParams,Functionals{1})
   plot(RegParams,Functionals{2})
end


end