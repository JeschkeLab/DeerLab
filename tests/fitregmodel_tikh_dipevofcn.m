function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.1;
RegMatrix = regoperator(Dimension,2);
TikhResult1 = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
TikhResult2 = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','bppnnls');
TikhResult3 = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','lsqnonneg','TolFun',1e-25);

err(1) = any(abs(TikhResult1 - P)>1e-4);
err(2) = any(abs(TikhResult2 - P)>1e-4);
err(3) = any(abs(TikhResult3 - P)>1e-4);

maxerr = max(abs(TikhResult1 - P));


err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,TikhResult1,'r')
end

end