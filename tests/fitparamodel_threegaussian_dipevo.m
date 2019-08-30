function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
TimeStep = 0.016;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [3 0.3 5 0.3 5.5 0.5 0.3 0.2];
Distribution = rd_threegaussian(DistanceAxis,InputParam);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

[FitDistribution] = fitparamodel(DipEvoFcn,Kernel,DistanceAxis,@rd_threegaussian,[],'solver','lsqnonlin');
err(1) = any(abs(FitDistribution - Distribution)>1e-1);
err = any(err);

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