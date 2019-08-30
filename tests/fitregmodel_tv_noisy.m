function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
NoiseLevel = 0.01;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Noise = whitenoise(Dimension,NoiseLevel);
Signal = DipEvoFcn + Noise;

RegMatrix = regoperator(Dimension,2);
RegParam = 0.001259;
TikhResult1 = fitregmodel(Signal,Kernel,DistanceAxis,RegMatrix,'tv',RegParam,'Solver','fmincon');
err(1) = any(abs(TikhResult1 - Distribution)>1e-1);
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