function [err,data] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(TimeAxis,DistanceAxis);

DipEvoFcn = Kernel*Distribution;
rng(2);
Noise = rand(Dimension,1);
Noise = Noise - mean(Noise);
Noise = 0.25*Noise;
NoiseLevel = std(Noise);
Signal = DipEvoFcn+Noise;

%Set optimal regularization parameter (found numerically lambda=0.13)
options = DAoptions('RegParam',500,'Solver','fmincon');
Result1 = OBIR(Signal,Kernel,'tikhonov',NoiseLevel,options);
err(1) = any(abs(Result1 - Distribution)>1e-2);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,Result1,'b')
end

end