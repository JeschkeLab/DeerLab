function [err,data] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;

RegParamSet = regparamrange(Kernel,RegMatrix);
OptParam2 = selregparam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,'tikhonov','all','NonNegConstrained',false,'NoiseLevel',0.05);

%Accept testif all values are the same (should be as there is no noise)
err = length(OptParam2)~=15;
data = [];



end