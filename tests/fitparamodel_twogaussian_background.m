function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [2.5 0.5 4 0.5 0.4];
Distribution = rd_twogaussian(DistanceAxis,InputParam);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Background = exp(-0.15*TimeAxis)';
ClusterFcn = (DipEvoFcn + 5).*Background;
Background = Background*(1-1/ClusterFcn(1));
ClusterFcn = ClusterFcn/ClusterFcn(1);
ClusterFcn = ClusterFcn./sqrt(Background);

KernelB = dipolarkernel(TimeAxis,DistanceAxis,Background,'KernelBType','sqrt');

InitialGuess = [2 0.1 5 0.1 0.1];
[FitDistribution] = fitparamodel(ClusterFcn,KernelB,DistanceAxis,@rd_twogaussian,InitialGuess);
err(1) = any(abs(FitDistribution - Distribution)>1e-5);
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