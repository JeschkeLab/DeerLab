function [err,data,maxerr] = test(opt,oldata)


Dimension = 500;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [20 0.49];
R0=0.198; % 1.98 Å per residue
SquareDist = 6*R0*InputParam(1)^InputParam(2)^2; %mean square end-to-end distance from radius of gyration
normFact = 3/(2*pi*SquareDist)^(3/2); % normalization prefactor
ShellSurf = 4*pi*r.^2; % spherical shell surface
Gaussian = exp(-3*r.^2/(2*SquareDist));
Distribution = normFact*ShellSurf.*Gaussian;
Distribution = Distribution/sum(Distribution)/mean(diff(r));
Distribution = Distribution.';

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

[FitDistribution] = fitparamodel(DipEvoFcn,K,r,@rd_randcoil);
err = any(abs(FitDistribution - Distribution)>1e-10);

maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf,
   subplot(121),hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitDistribution,'r')
   subplot(122),hold on
   plot(r,Distribution,'b')
   plot(r,FitDistribution,'r')
end

end