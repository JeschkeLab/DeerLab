function [err,data,maxerr] = test(opt,olddata)

Ntime = 100;
Ndist = 100;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);

Distribution = rd_onegaussian(r,[3,0.5]);
Distribution = Distribution/sum(Distribution);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;
RegMatrix = regoperator(Ndist,2);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParamRange = regparamrange(K,RegMatrix);
RegParam1 = selregparam(RegParamRange,DipEvoFcn,K,RegMatrix,'tikhonov','aic');


Ntime = 100;
Ndist = 200;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);

Distribution = rd_onegaussian(r,[3,0.5]);
Distribution = Distribution/sum(Distribution);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;
RegMatrix = regoperator(Ndist,2);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParamRange = regparamrange(K,RegMatrix);
RegParam2 = selregparam(RegParamRange,DipEvoFcn,K,RegMatrix,'tikhonov','aic');

%RegParam2 should be larger to compensate for worse condition number 
err = RegParam2 < RegParam1;

maxerr = [];
data = [];

end