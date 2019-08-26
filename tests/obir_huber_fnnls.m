function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = twogaussian(DistanceAxis,[2,0.3,3.5,0.3,0.5]);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix =  regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;
NoiseLevel = 0.05;
Noise = whitenoise(Dimension,NoiseLevel);
Signal = DipEvoFcn+Noise;

%Set optimal regularization parameter (found numerically lambda=0.13)
OptParam = 12;
OptHuber = 1.35;

if opt.Display
    figure(8),clf
    axhandle = plot(DistanceAxis,NaN*Distribution);
else
    axhandle = [];
end

Result = obir(Signal,Kernel,DistanceAxis,'huber',RegMatrix,OptParam,'DivergenceStop',true,'NoiseLevelAim',NoiseLevel,'Solver','fnnls','Huberparam',OptHuber,'axishandle',axhandle);

RegResult = fitregmodel(Signal,Kernel,DistanceAxis,RegMatrix,'huber',OptParam);

err = norm(Result - Distribution) > norm(RegResult - Distribution);
maxerr = norm(Result - Distribution);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,Result,'b')
    plot(DistanceAxis,RegResult,'r')
    legend('truth','OBIR','Huber')
end

end