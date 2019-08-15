function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 1e-3;
RegMatrix = regoperator(Dimension,3);
TikhResult1 = fitregmodel(DipEvoFcn,DistanceAxis,Kernel,RegMatrix,'tv',RegParam,'Solver','fnnls');

err = any(abs(TikhResult1 - Distribution)>2e-2);

maxerr = max(abs(TikhResult1 - Distribution));

data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TikhResult1,'r')
    axis tight
end

end