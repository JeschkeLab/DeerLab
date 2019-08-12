function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [2 0.5 3 0.5 0.4];
Distribution = InputParam(5)*gaussian(DistanceAxis,InputParam(1),InputParam(2)) + (1 - InputParam(5))*gaussian(DistanceAxis,InputParam(3),InputParam(4));
Distribution = Distribution/(1/sqrt(2*pi)*1/InputParam(2));
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

InitialGuess = [2 0.1 5 0.1 0.5];
[FitDistribution,FitParam] = fitparamodel(DipEvoFcn,Kernel,DistanceAxis,@twogaussian,InitialGuess);
[multigaussFitDistribution,~,N] = multigauss(DipEvoFcn,Kernel,DistanceAxis,2);

err(1) = any(abs(FitDistribution - multigaussFitDistribution)>1e-9);
err = any(err);

maxerr = max(abs(FitDistribution - multigaussFitDistribution));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(TimeAxis,Kernel*FitDistribution,'b')
   plot(TimeAxis,Kernel*multigaussFitDistribution,'r')
   subplot(122)
   hold on
   plot(DistanceAxis,FitDistribution,'b')
   plot(DistanceAxis,multigaussFitDistribution,'r')
   
end

end