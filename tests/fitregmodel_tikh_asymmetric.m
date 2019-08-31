function [err,data,maxerr] = test(opt,olddata)

Ntime = 100;
Ndist = 200;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;
RegMatrix = regoperator(Ndist,2);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParamRange = regparamrange(K,RegMatrix);
RegParam = selregparam(RegParamRange,DipEvoFcn,K,RegMatrix,'tikhonov','aic');

TikhResult = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');


err(1) = any(abs(TikhResult - Distribution)>3e-3);
err(2) = length(TikhResult) ~= Ndist;
err(3) = length(K*TikhResult) ~= Ntime;
err  = any(err);

maxerr = max(abs(TikhResult - Distribution));


err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,Distribution,'k') 
    plot(r,TikhResult,'r')
end

end