function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
P1 = onegaussian(DistanceAxis,[2,0.3]);
P2 = onegaussian(DistanceAxis,[3.5,0.3]);
Distribution = 0.5*P1 + 0.5*P2;
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix =  regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;
NoiseLevel = 0.05;
Noise = whitenoise(Dimension,NoiseLevel);
Signal = DipEvoFcn+Noise;

if opt.Display
    figure(8),clf
    axhandle = plot(DistanceAxis,NaN*Distribution);
else
    axhandle = [];
end

%Set optimal regularization parameter (found numerically lambda=0.13)
OptParam = 0.1;
Result = obir(Signal,Kernel,DistanceAxis,'tv',RegMatrix,OptParam,'DivergenceStop',true,'NoiseLevelAim',NoiseLevel,'Solver','fnnls','axishandle',axhandle);
RegResult = fitregmodel(Signal,Kernel,DistanceAxis,RegMatrix,'tv',OptParam);

err = norm(Result - Distribution) > norm(RegResult - Distribution);
maxerr = norm(Result - Distribution);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,Result,'b')
    plot(DistanceAxis,RegResult,'r')
    legend('truth','OBIR','TV')
end

end