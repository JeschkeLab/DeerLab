function [err,data,maxerr] = test(opt,olddata)

Ntime1 = 100;
Ndist = 200;

TimeStep = 0.008;
TimeAxis1 = linspace(0,TimeStep*Ntime1,Ntime1);
[~,rmin,rmax] = time2dist(TimeAxis1);
DistanceAxis = linspace(rmin,rmax,Ndist);

Distribution = twogaussian(DistanceAxis,[2,0.3,4,0.3,0.5]);

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

Result = fitregmodel(Signals,Kernels,DistanceAxis,L,'tikhonov',regparam,'Solver','fnnls');
Dist1 = fitregmodel(Signal1,Kernel1,DistanceAxis,L,'tikhonov',regparam,'Solver','fnnls');
Dist2 = fitregmodel(Signal2,Kernel2,DistanceAxis,L,'tikhonov',regparam,'Solver','fnnls');
Dist3 = fitregmodel(Signal3,Kernel3,DistanceAxis,L,'tikhonov',regparam,'Solver','fnnls');

normResult = norm(Distribution - Result);
norm1 = norm(Distribution - Dist1);
norm2 = norm(Distribution - Dist2);
norm3 = norm(Distribution - Dist3);

err = any(normResult > [norm2 norm3]);
err = any(err);
data = [];
maxerr = normResult;

if opt.Display
figure(8),clf
subplot(121)
hold on
plot(TimeAxis1,Signal1,'k')
plot(TimeAxis2,Signal2+1,'k')
plot(TimeAxis3,Signal3+2,'k')
plot(TimeAxis1,Kernel1*Result,'r')
plot(TimeAxis2,Kernel2*Result + 1,'r')
plot(TimeAxis3,Kernel3*Result + 2,'r')
subplot(122)
hold on
plot(DistanceAxis,Distribution,'k')
plot(DistanceAxis,Result,'r')
plot(DistanceAxis,Dist1,'g--')
plot(DistanceAxis,Dist2,'b--')
plot(DistanceAxis,Dist3,'r--')
end

end