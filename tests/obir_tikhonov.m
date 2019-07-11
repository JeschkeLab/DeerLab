function [err,data] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
TimeStep = 0.008;
rmin = (4*TimeStep*52.04/0.85)^(1/3);
rmax = 6*(Dimension*TimeStep/2)^(1/3);
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = linspace(rmin,rmax,Dimension);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(Dimension,TimeStep*1000);
Background = exp(-0.5*TimeAxis);

DipEvoFcn = Kernel*Distribution;
rng(2);
Noise = rand(Dimension,1);
Noise = Noise - mean(Noise);
Noise = 0.25*Noise;
NoiseLevel = std(Noise);
Signal = DipEvoFcn+Noise;

%Set optimal regularization parameter (found numerically lambda=0.13)
options = DAoptions('RegParam',1000,'Solver','fmincon');
Result1 = OBIR(Signal,Kernel,'tikhonov',NoiseLevel,options);

Result2 = regularize(Signal,Kernel,'tikhonov',options);
err(1) = any(abs(Result1 - Distribution)>1e-2);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,Result1,'b')
    plot(DistanceAxis,Result2,'r')

end

end