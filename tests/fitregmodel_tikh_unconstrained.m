function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 1;
Result = fitregmodel(DipEvoFcn,K,r,'tikhonov',RegParam,'NonNegConstrained',false);

err(1) = any(abs(Result - P)>4e-2);
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