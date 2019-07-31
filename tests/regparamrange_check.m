function [err,data] = test(opt,olddata)

%=======================================
% Check regparamrange.m
%=======================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
L = regoperator(Dimension,2);
logRes = 0.1;
alpha = regparamrange(Kernel,L,'logResolution',logRes);

err = (length(alpha)~=85);
data = [];



end