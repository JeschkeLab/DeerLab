function [pass,maxerr] = test(opt)

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
RegParam = 1e-3;
TikhResult1 = fitregmodel(DipEvoFcn,K,r,'tv',RegParam,'Solver','fnnls','RegOrder',3);

pass = all(abs(TikhResult1 - P)>3e-2);

maxerr = max(abs(TikhResult1 - P));

 

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,TikhResult1,'r')
    axis tight
end

end