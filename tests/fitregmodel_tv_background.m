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
V = (DipEvoFcn + 5).*B;
ModDepth = 1/V(1);
V = V/V(1);


%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.0005;
KB = dipolarkernel(t,r,ModDepth,B);
Result = fitregmodel(V,KB,r,'tv',RegParam,'Solver','fnnls','RegOrder',3);

error = abs(Result - P);
err(1) = any(error > 2e-2);
maxerr = max(error);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'r')
end

end