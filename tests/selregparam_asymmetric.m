function [err,data,maxerr] = test(opt,olddata)

Ntime = 100;
Ndist = 100;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);

P = rd_onegaussian(r,[3,0.5]);
P = P/sum(P);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam1 = selregparam(DipEvoFcn,K,'tikhonov','aic');


Ntime = 100;
Ndist = 200;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);

P = rd_onegaussian(r,[3,0.5]);
P = P/sum(P);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam2 = selregparam(DipEvoFcn,K,'tikhonov','aic');

%RegParam2 should be larger to compensate for worse condition number 
err = RegParam2 < RegParam1;

maxerr = [];
data = [];

end