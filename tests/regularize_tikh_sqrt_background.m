function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussfcn(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Background = exp(-0.15*TimeAxis)';
ClusterFcn = (DipEvoFcn + 5).*Background;
Background = Background*(1-1/ClusterFcn(1));
ClusterFcn = ClusterFcn/ClusterFcn(1);
ClusterFcn = ClusterFcn./sqrt(Background);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.13;
KernelB = dipolarkernel(TimeAxis,DistanceAxis,Background,'KernelBType','sqrt');
RegMatrix = regoperator(Dimension,2);
TikhResult1 = regularize(ClusterFcn,KernelB,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
TikhResult2 = regularize(ClusterFcn,KernelB,RegMatrix,'tikhonov',RegParam,'Solver','bppnnls');
TikhResult3 = regularize(ClusterFcn,KernelB,RegMatrix,'tikhonov',RegParam,'Solver','lsqnonneg','nonNegLSQsolTol',1e-29);

err(1) = any(abs(TikhResult1 - Distribution)>1e-5);
err(2) = any(abs(TikhResult2 - Distribution)>1e-5);
err(3) = any(abs(TikhResult3 - Distribution)>1e-2);

maxerr = max(abs(TikhResult1 - Distribution));
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TikhResult1,'r')
end

end