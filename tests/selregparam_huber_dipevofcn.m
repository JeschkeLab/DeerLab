function [err,data] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 80;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = onegaussian(DistanceAxis,[3,0.5]);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;

RegParamSet = regparamrange(Kernel,RegMatrix);
[OptParam,~,~] = selregparam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,'huber',{'aic','gcv'});

%Accept testif all values are the same (should be as there is no noise)
err = any(any(OptParam - OptParam' > 1e-2));
data = [];



end