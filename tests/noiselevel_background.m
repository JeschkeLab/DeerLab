function [err,data,maxerr] = test(opt,olddata)

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [3 0.5];
Distribution = onegaussian(DistanceAxis,InputParam);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;
Background = exp(-1.5*TimeAxis)';
ClusterFcn = (DipEvoFcn + 5).*Background;
Background = Background*(1-1/ClusterFcn(1));
ClusterFcn = ClusterFcn/ClusterFcn(1);
ClusterFcn = ClusterFcn./sqrt(Background);

rng(2)
Noise = rand(Dimension,1);
Noise = Noise - mean(Noise);
Noise = 0.02*Noise/Noise(1);

ClusterFcn = ClusterFcn + Noise;

truelevel = std(Noise);
approxlevel = noiselevel(ClusterFcn);

err = abs(approxlevel - truelevel)>1e-2;
maxerr = abs(approxlevel - truelevel);
data = [];


end
