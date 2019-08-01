function [err,data,maxerr] = test(opt,oldata)


Dimension = 300;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [2 0.4 3.5 0.3 5 0.3 0.3 0.3];
Distribution = threerice(DistanceAxis,InputParam);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

[FitDistribution] = fitparamodel(DipEvoFcn,Kernel,DistanceAxis,@threerice,[],'Constrained',true);
err = any(abs(FitDistribution - Distribution)>1e-5);

maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(TimeAxis,DipEvoFcn,'b')
   plot(TimeAxis,Kernel*FitDistribution,'r')
   subplot(122)
   hold on
   plot(DistanceAxis,Distribution,'b')
   plot(DistanceAxis,FitDistribution,'r')
end

end