function [err,data] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 80;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;

RegParamSet = regparamrange(Kernel,RegMatrix);
OptParam = selregparam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,{'aic','gcv','lr'});

%Accept testif all values are the same (should be as there is no noise)
err(1) = any(diff(OptParam) > 1e-2);
err(2) = any(abs(OptParam - 0.001995262314969));
err = any(err);
data = [];



end