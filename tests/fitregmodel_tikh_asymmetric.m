function [err,data,maxerr] = test(opt,olddata)

Ntime = 100;
Ndist = 200;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

TikhResult = fitregmodel(DipEvoFcn,K,r,'tikhonov','aic','Solver','fnnls');


err(1) = any(abs(TikhResult - P)>3e-2);
err(2) = length(TikhResult) ~= Ndist;
err(3) = length(K*TikhResult) ~= Ntime;
err  = any(err);

maxerr = max(abs(TikhResult - P));


err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,TikhResult,'r')
end

end