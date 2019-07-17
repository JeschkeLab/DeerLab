function [err,data] = test(opt,olddata)

%=======================================
% Check getRegParamRange.m
%=======================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2distAxis(TimeAxis);

Kernel = getKernel(TimeAxis,DistanceAxis);
L = getRegMatrix(Dimension,2);
logRes = 0.1;
alpha = getRegParamRange(Kernel,L,'logResolution',logRes);

err = (length(alpha)~=85);
data = [];



end