function [err,data,maxerr] = test(opt,oldata)


Dimension = 500;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [20 0.49];
R0=0.198; % 1.98 Å per residue
SquareDist = 6*R0*InputParam(1)^InputParam(2)^2; %mean square end-to-end distance from radius of gyration
normFact = 3/(2*pi*SquareDist)^(3/2); % normalization prefactor
ShellSurf = 4*pi*DistanceAxis.^2; % spherical shell surface
Gaussian = exp(-3*DistanceAxis.^2/(2*SquareDist));
Distribution = normFact*ShellSurf.*Gaussian;
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));
Distribution = Distribution.';

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

[FitDistribution] = fitparamodel(DipEvoFcn,Kernel,DistanceAxis,@randcoil);
err = any(abs(FitDistribution - Distribution)>1e-10);

maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf,
   subplot(121),hold on
   plot(TimeAxis,DipEvoFcn,'b')
   plot(TimeAxis,Kernel*FitDistribution,'r')
   subplot(122),hold on
   plot(DistanceAxis,Distribution,'b')
   plot(DistanceAxis,FitDistribution,'r')
end

end