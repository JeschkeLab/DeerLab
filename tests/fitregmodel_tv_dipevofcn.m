function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check TV regularization
%=======================================
Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

%Set optimal regularization parameter (found numerically lambda=0.005)
RegParam = 0.0005;
RegMatrix = regoperator(Dimension,3);
TVResult1 = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tv',RegParam,'Solver','fmincon');

error = abs(TVResult1 - P);
err(1) = any(error>2e-1);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,TVResult1,'r')
end

end