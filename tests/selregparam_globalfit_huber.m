function [err,data,maxerr] = test(opt,olddata)

Ntime1 = 100;
Ndist = 200;
TimeStep = 0.008;
TimeAxis1 = linspace(0,TimeStep*Ntime1,Ntime1);
[~,rmin,rmax] = time2dist(TimeAxis1);
DistanceAxis = linspace(rmin,rmax,Ndist);

Distribution = gaussian(DistanceAxis,2,0.3) + gaussian(DistanceAxis,3.5,0.3);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

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

Signals = {Signal1,Signal2,Signal3};
Kernels = {Kernel1,Kernel2,Kernel3};


RegParamRange = logspace(-3,4,60);
[OptRegParam,fun] = selregparam(RegParamRange,Signals,Kernels,L,'huber',{'aic','aicc'});

Result1 = fitregmodel(Signals,Kernels,DistanceAxis,L,'huber',OptRegParam(1),'Solver','fnnls');
Result2 = fitregmodel(Signals,Kernels,DistanceAxis,L,'huber',OptRegParam(2),'Solver','fnnls');

Result = mean([Result1 Result2],2);
stdDist = std([Result1 Result2],1,2);

err = any(abs(Result-Distribution) > 0.5);
data = [];
maxerr = max(abs(Result-Distribution));

if opt.Display
figure(8),clf
subplot(121)
hold on
plot(RegParamRange,fun{1},'.b')
subplot(122)
hold on
plot(DistanceAxis,Distribution,'k')
Result = Result';
stdDist = stdDist';
f = fill([DistanceAxis fliplr(DistanceAxis)],[Result+stdDist fliplr(Result-stdDist)],'b','LineStyle','none');
f.FaceAlpha = 0.5;
plot(DistanceAxis,Result,'b')
end

end