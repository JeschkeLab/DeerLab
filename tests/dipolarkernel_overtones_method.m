
function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check overtones in kernel calculation
%======================================================

N = 50;
t = linspace(0,3,N); % us
r = time2dist(t); % nm
Tmix = 50; % us
T1 = 88; % us

coefficients = overtones(3,Tmix,T1);

Kfresnel = dipolarkernel(t,r,'Method','fresnel','OvertoneCoeffs',coefficients);
Kgrid = dipolarkernel(t,r,'Method','grid','OvertoneCoeffs',coefficients);

dr = mean(diff(r));
delta = abs(Kfresnel-Kgrid)/dr;
err = any(delta(:)>1e-3);
maxerr = max(delta(:));

data = [];

end