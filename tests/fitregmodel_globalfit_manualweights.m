function [err,data,maxerr] = test(opt,olddata)

Ntime1 = 100;
Ndist = 200;

TimeStep = 0.008;
TimeAxis1 = linspace(0,TimeStep*Ntime1,Ntime1);
[~,rmin,rmax] = time2dist(TimeAxis1);
DistanceAxis = linspace(rmin,rmax,Ndist);

Distribution = rd_twogaussian(DistanceAxis,[2,0.3,4,0.3,0.5]);

Kernel1 = dipolarkernel(TimeAxis1,DistanceAxis);
Signal1 = Kernel1*Distribution;
noise = whitenoise(length(Signal1),0.03);
Signal1 = Signal1 + noise;

Ntime2 = 200;
TimeAxis2 = linspace(0,TimeStep*Ntime2,Ntime2);
Kernel2 = dipolarkernel(TimeAxis2,DistanceAxis);
Signal2 = Kernel2*Distribution;
noise = whitenoise(length(Signal2),0.05);
Signal2 = Signal2 + noise;

Ntime3 = 300;
TimeAxis3 = linspace(0,TimeStep*Ntime3,Ntime3);
Kernel3 = dipolarkernel(TimeAxis3,DistanceAxis);
Signal3 = Kernel3*Distribution;
noise = whitenoise(length(Signal3),0.1);
Signal3 = Signal3 + noise;


L = regoperator(Ndist,2);
%Set optimal regularization parameter (found numerically lambda=0.13)
regparam = 2;

Signals = {Signal1,Signal2,Signal3};
Kernels = {Kernel1,Kernel2,Kernel3};

Result = fitregmodel(Signals,Kernels,DistanceAxis,L,'tikhonov',regparam,'Solver','fnnls','GlobalWeights',[0.4 0.5 0.1]);

err = any(isnan(Result));
err = any(err);
data = [];
maxerr = [];

end