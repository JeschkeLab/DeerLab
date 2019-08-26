function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%gaussian model as given as the input signal

Dimension = 500;
TimeStep = 0.016;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [3 0.3 5 0.3 0.5];
Distribution = twogaussian(DistanceAxis,InputParam);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

Models = {@onegaussian,@twogaussian,@threegaussian,@onerice,@tworice};

tic
[preoptimum] = selectmodel(Models,DipEvoFcn,DistanceAxis,Kernel,'aicc');
pre  = toc;
tic
[postoptimum] = selectmodel(Models,DipEvoFcn,DistanceAxis,Kernel,'aicc');
post  = toc;

err(1) = any(postoptimum~=preoptimum);
err(2) = post > pre/3;
data = [];
maxerr = [];


end

