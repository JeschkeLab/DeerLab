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
RegParam = 5;
Result = fitregmodel(DipEvoFcn,K,r,'tikhonov',RegParam,'Solver','fmincon');

err(1) = any(abs(Result - P)>5e-1);
maxerr = max(abs(Result - P));
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'r')
end

end