function [err,data,maxerr] = test(opt,data)

%Clear persistent variable
clear regularize

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;

RegParam = 100;

tic
preDist = fitregmodel(DipEvoFcn,Kernel,DistanceAxis,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
precached = toc;

tic
postDist = fitregmodel(DipEvoFcn,Kernel,DistanceAxis,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
postcached = toc;

err(1) = postcached>=precached/4;
error = abs(preDist - postDist);
err(2) = any(any(error>1e-18));
err = any(err);
data = [];
maxerr = max(max(error));

end