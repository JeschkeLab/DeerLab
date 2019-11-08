function [err,data,maxerr] = test(opt,olddata)

clear regularize

%=======================================
% Check TV regularization
%=======================================
Dimension = 80;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

%Set optimal regularization parameter (found numerically lambda=0.005)
RegParam = 1;
Result = fitregmodel(DipEvoFcn,K,r,'huber',RegParam,'Solver','fmincon','HuberParam',1.35);

error = abs(Result - P);
err(1) = any(error>2e-1);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'r')
end

end