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
B = exp(-0.15*t)';
V = (DipEvoFcn+5).*B;
ModDepth = 1/V(1);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.13;
KB = dipolarkernel(t,r,ModDepth,B);
TikhResult1 = fitregmodel(V,KB,r,'tikhonov',RegParam,'Solver','fnnls');
TikhResult2 = fitregmodel(V,KB,r,'tikhonov',RegParam,'Solver','bppnnls');
TikhResult3 = fitregmodel(V,KB,r,'tikhonov',RegParam,'Solver','lsqnonneg','TolFun',1e-15);

err(1) = any(abs(TikhResult1 - P)>3e-3);
err(2) = any(abs(TikhResult2 - P)>3e-3);
err(3) = any(abs(TikhResult3 - P)>3e-3);

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