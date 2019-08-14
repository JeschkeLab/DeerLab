function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
NoiseLevel = 0.01;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Noise = whitenoise(Dimension,NoiseLevel);
Signal = DipEvoFcn + Noise;

RegMatrix = regoperator(Dimension,2);
range = regparamrange(Kernel,RegMatrix);
[RegParam] = selregparam(range,Signal,Kernel,RegMatrix,'aic','RegType','huber');

TikhResult1 = regularize(Signal,DistanceAxis,Kernel,RegMatrix,'huber',RegParam,'Solver','fnnls');
err(1) = any(abs(TikhResult1 - Distribution)>5e-2);
maxerr = max(abs(TikhResult1 - Distribution));

err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    subplot(121)
    hold on
    plot(TimeAxis,Signal)
    plot(TimeAxis,Kernel*TikhResult1)
    subplot(122)
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TikhResult1,'r')
end

end