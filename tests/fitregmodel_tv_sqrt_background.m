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
V = V./sqrt(B);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.005;
RegMatrix = regoperator(Dimension,3);
KB = dipolarkernel(t,r,ModDepth,sqrt(B));
Result = fitregmodel(V,KB,r,RegMatrix,'tv',RegParam,'Solver','fnnls');

err = any(abs(Result - P)>1.5e-2);
maxerr = max(abs(Result - P));

data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'r')
end

end