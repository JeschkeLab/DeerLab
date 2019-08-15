function [err,data,maxerr] = test(opt,olddata)

Ntime = 100;
Ndist = 200;

TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Ntime,Ntime);
[~,rmin,rmax] = time2dist(TimeAxis);
DistanceAxis = linspace(rmin,rmax,Ndist);

Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
RegMatrix = regoperator(Ndist,2);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParamRange = regparamrange(Kernel,RegMatrix);
RegParam = selregparam(RegParamRange,DipEvoFcn,Kernel,RegMatrix,'aic');

TikhResult = fitregmodel(DipEvoFcn,DistanceAxis,Kernel,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');


err(1) = any(abs(TikhResult - Distribution)>3e-3);
err(2) = length(TikhResult) ~= Ndist;
err(3) = length(Kernel*TikhResult) ~= Ntime;
err  = any(err);

maxerr = max(abs(TikhResult - Distribution));


err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TikhResult,'r')
end

end