function [err,data] = test(opt,olddata)

%=======================================
% Check getRegParamRange.m
%=======================================

Dimension = 80;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(TimeAxis,DistanceAxis);
RegMatrix = getRegMatrix(Dimension,2);
DipEvoFcn = Kernel*Distribution;

RegParamSet = getRegParamRange(Kernel,RegMatrix);
OptParam = selectRegParam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,{'aic','gcv'},'RegType','tv');

%Accept testif all values are the same (should be as there is no noise)
err = any(any(OptParam - OptParam' > 1e-2));
data = [];



end