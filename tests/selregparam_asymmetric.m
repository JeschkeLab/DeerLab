function [err,data,maxerr] = test(opt,olddata)

Ntime = 100;
Ndist = 100;

TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Ntime,Ntime);
[~,rmin,rmax] = time2dist(TimeAxis);
DistanceAxis = linspace(rmin,rmax,Ndist);

Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
RegMatrix = regoperator(Ndist,2);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParamRange = regparamrange(Kernel,RegMatrix);
RegParam1 = selregparam(RegParamRange,DipEvoFcn,Kernel,RegMatrix,'aic');


Ntime = 100;
Ndist = 200;

TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Ntime,Ntime);
[~,rmin,rmax] = time2dist(TimeAxis);
DistanceAxis = linspace(rmin,rmax,Ndist);

Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
RegMatrix = regoperator(Ndist,2);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParamRange = regparamrange(Kernel,RegMatrix);
RegParam2 = selregparam(RegParamRange,DipEvoFcn,Kernel,RegMatrix,'aic');

%RegParam2 should be larger to compensate for worse condition number 
err = RegParam2 < RegParam1;

maxerr = [];
data = [];

end