function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [2 0.3];
Distribution = onegaussian(DistanceAxis,[2,0.3]);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

[FitDistribution] = fitparamodel(DipEvoFcn,Kernel,DistanceAxis,@onegaussian,0.6*InputParam,'Solver','fmincon','costmodel','chisquare');
err(1) = any(abs(FitDistribution - Distribution)>1e-7);
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
   legend('truth','fit')
end

end