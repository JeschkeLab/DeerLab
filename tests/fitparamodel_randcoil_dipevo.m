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
P = normFact*ShellSurf.*Gaussian;
P = P/sum(P)/mean(diff(r));
P = P.';

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

[~,FitP] = fitparamodel(DipEvoFcn,@rd_randcoil,r,K);
err = any(abs(FitP - P)>1e-10);

maxerr = max(abs(FitP - P));
data = [];

if opt.Display
   figure(1),clf,
   subplot(121),hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitP,'r')
   subplot(122),hold on
   plot(r,P,'b')
   plot(r,FitP,'r')
end

end