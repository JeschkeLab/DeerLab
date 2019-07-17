function [err,data] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.1;
TikhResult1 = regularize(DipEvoFcn,Kernel,'tikhonov',RegParam,'Solver','fnnls');
TikhResult2 = regularize(DipEvoFcn,Kernel,'tikhonov',RegParam,'Solver','bppnnls');
TikhResult3 = regularize(DipEvoFcn,Kernel,'tikhonov',RegParam,'Solver','lsqnonneg','nonNegLSQsolTol',1e-25);

err(1) = any(abs(TikhResult1 - Distribution)>1e-6);
err(2) = any(abs(TikhResult2 - Distribution)>1e-6);
err(3) = any(abs(TikhResult3 - Distribution)>1e-5);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TikhResult1,'r')
end

end