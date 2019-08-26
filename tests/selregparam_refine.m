function [err,data,maxerr] = test(opt,olddata)
Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;


Noise = rand(Dimension,1);
Noise = Noise - mean(Noise);
Noise = 0.05*Noise/max(Noise);

RegParamSet = regparamrange(Kernel,RegMatrix);
goal = 0.01;
[~,Functionals,RegParams] = selregparam(RegParamSet,DipEvoFcn + Noise,Kernel,RegMatrix,'tikhonov',{'aic','gcv'},'Refine',true);

err = length(RegParams) == length(RegParamSet);
data = [];
maxerr = [];

if opt.Display
   figure(8),clf
hold on
a = Functionals{1};
plot((RegParams(1:end-60)),a(1:end-60),'.')
plot((RegParams(end-60:end)),a(end-60:end),'.')
end

end