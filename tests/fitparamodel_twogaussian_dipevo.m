function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [2 0.5 3 0.5 0.4];
Distribution = rd_twogaussian(DistanceAxis,InputParam);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

InitialGuess = [2 0.1 5 0.1 0.5];
[FitDistribution,FitParam] = fitparamodel(DipEvoFcn,Kernel,DistanceAxis,@rd_twogaussian,InitialGuess);
err(1) = any(abs(FitDistribution - Distribution)>1e-5);
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