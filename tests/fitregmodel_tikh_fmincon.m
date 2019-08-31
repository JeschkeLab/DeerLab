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
RegParam = 5;
RegMatrix = regoperator(Dimension,2);
Result = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fmincon');

err(1) = any(abs(Result - Distribution)>5e-1);
maxerr = max(abs(Result - Distribution));
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,Distribution,'k') 
    plot(r,Result,'r')
end

end