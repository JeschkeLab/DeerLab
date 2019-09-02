function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.1;
RegMatrix = regoperator(Dimension,2);
%nm
Pfit1 = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
%A
r = r*10;
Pfit2 = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');

err = any(abs(Pfit1 - Pfit2)>1e-12);
maxerr = max(abs(Pfit1 - Pfit2));
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,Pfit1,'k') 
    plot(r,Pfit2,'r')
end

end