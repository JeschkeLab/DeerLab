function [err,data,maxerr] = test(opt,oldata)

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [3 0.5];
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Background = exp(-0.15*TimeAxis)';
ClusterFcn = (DipEvoFcn + 5).*Background;
Background = Background*(1-1/ClusterFcn(1));
ClusterFcn = ClusterFcn/ClusterFcn(1);
ClusterFcn = ClusterFcn./sqrt(Background);

KernelB = dipolarkernel(TimeAxis,DistanceAxis,Background,'KernelBType','sqrt');

InitialGuess = [2 0.1];
[FitDistribution,FitParam] = fitparamodel(ClusterFcn,KernelB,DistanceAxis,@rd_onegaussian,InitialGuess);
err(1) = any(abs(FitDistribution - Distribution)>1e-5);
err(2) = any(abs(FitParam - InputParam)>1e-3);
err = any(err);

maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(TimeAxis,ClusterFcn,'b')
   plot(TimeAxis,KernelB*FitDistribution,'r')
end

end