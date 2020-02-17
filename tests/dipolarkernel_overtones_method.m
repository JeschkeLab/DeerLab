
function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check overtones in kernel calculation
%======================================================

N = 100;
t = linspace(0,3,N);
r = time2dist(t);
Tmix = 50; % us
T1 = 88; % us

coefficients = overtones(3,Tmix,T1);

Kfresnel = dipolarkernel(t,r,'Method','fresnel','OvertoneCoeffs',coefficients);
Kgrid = dipolarkernel(t,r,'Method','grid','OvertoneCoeffs',coefficients);

dr = mean(diff(r));
Kfresnel = Kfresnel/dr;
Kgrid = Kgrid/dr;

delta = abs(Kfresnel - Kgrid);
err = any(delta(:)>1e-3);
maxerr = max(delta(:));

data = [];

end