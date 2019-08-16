function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Noise = whitenoise(Dimension,0.02);
%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.1;
RegMatrix = regoperator(Dimension,3);
Resultfnnls = fitregmodel(DipEvoFcn+Noise,Kernel,DistanceAxis,RegMatrix,'tv',RegParam,'Solver','fnnls');

err = any(abs(Resultfnnls - Distribution)>9e-2);

maxerr = max(abs(Resultfnnls - Distribution));

data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k')
    plot(DistanceAxis,Resultfnnls,'b') 
    axis tight
end

end