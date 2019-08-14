function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [2.5 0.5 4 0.5 3 0.2 0.3 0.4];
Distribution = InputParam(7)*gaussian(DistanceAxis,InputParam(1),InputParam(2))/(1/sqrt(2*pi)*1/InputParam(2)) ...
    + InputParam(8)*gaussian(DistanceAxis,InputParam(3),InputParam(4))/(1/sqrt(2*pi)*1/InputParam(4));
    + (1 - InputParam(7) -InputParam(8))*gaussian(DistanceAxis,InputParam(5),InputParam(6))/(1/sqrt(2*pi)*1/InputParam(6));
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Background = exp(-0.15*TimeAxis)';
ClusterFcn = (DipEvoFcn + 5).*Background;
Background = Background*(1-1/ClusterFcn(1));
ClusterFcn = ClusterFcn/ClusterFcn(1);
ClusterFcn = ClusterFcn./sqrt(Background);

KernelB = dipolarkernel(TimeAxis,DistanceAxis,Background,'KernelBType','sqrt');

InitialGuess = [2 0.1 5 0.1 1 0.2 0.1 0.5];
[FitDistribution] = fitparamodel(ClusterFcn,KernelB,DistanceAxis,@threegaussian,InitialGuess);
err(1) = any(abs(FitDistribution - Distribution)>1e-4);
err = any(err);
maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(TimeAxis,ClusterFcn,'b')
   plot(TimeAxis,KernelB*FitDistribution,'r')
   subplot(122)
   hold on
   plot(DistanceAxis,Distribution,'b')
   plot(DistanceAxis,FitDistribution,'r')
   
end

end