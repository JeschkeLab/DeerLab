function [err,data,maxerr] = test(opt,oldata)

% Build time and distance axes
nt = 500;
dt = 0.008;
t = linspace(0,dt*nt,nt);
r = time2dist(t);

% Build random-coil model distribution
N = 20; % number of residues
nu = 0.49; % scaling exponent
InputParam = [N nu];
R0 = 0.198; % 1.98 angstrom per residue
rsq = 6*(R0*N^nu)^2; % mean square end-to-end distance from radius of gyration
normFact = 3/(2*pi*rsq)^(3/2); % normalization prefactor
ShellSurf = 4*pi*r.^2; % spherical shell surface
Gaussian = exp(-3*r.^2/(2*rsq));
P = normFact*ShellSurf.*Gaussian;
P = P/sum(P)/mean(diff(r)); % normalize
P = P(:);

% Calculate model dipolar signal
K = dipolarkernel(t,r);
V = K*P;

% Fit random-coil model to dipolar signal
[~,Pfit] = fitparamodel(V,@rd_randcoil,r,K);

err = any(abs(Pfit - P)>1e-7);
maxerr = max(abs(Pfit - P));
data = [];

if opt.Display
    figure(1),clf,
    subplot(121),hold on
    plot(t,DipEvoFcn,'b')
    plot(t,K*Pfit,'r')
    subplot(122),hold on
    plot(r,P,'b')
    plot(r,Pfit,'r')
end

end
