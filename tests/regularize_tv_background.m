function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Background = exp(-0.15*TimeAxis)';
ClusterFcn = (DipEvoFcn + 5).*Background;
Background = Background*(1-1/ClusterFcn(1));
ClusterFcn = ClusterFcn/ClusterFcn(1);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.005;
RegMatrix = regoperator(Dimension,3);
KernelB = dipolarkernel(TimeAxis,DistanceAxis,Background,'KernelBType','full');
Result = regularize(ClusterFcn,KernelB,RegMatrix,'tv',RegParam,'Solver','fmincon');

error = abs(Result - Distribution);
err(1) = any(error > 5e-3);
maxerr = max(error);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,Result,'r')
end

end