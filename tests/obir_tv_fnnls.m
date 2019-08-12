function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = 0.5*gaussian(DistanceAxis,2,0.3) + 0.5*gaussian(DistanceAxis,3.5,0.3);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix =  regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;
NoiseLevel = 0.05;
Noise = whitenoise(Dimension,NoiseLevel);
Signal = DipEvoFcn+Noise;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParamSet = regparamrange(Kernel,RegMatrix);
[OptParam] = selregparam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,'gcv','RegType','tv');
Result = obir(Signal,Kernel,'tv',RegMatrix,OptParam,'DivergenceStop',true,'NoiseLevelAim',NoiseLevel,'Solver','fnnls');

RegResult = regularize(Signal,Kernel,RegMatrix,'tv',OptParam);

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