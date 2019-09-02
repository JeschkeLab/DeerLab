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
B = exp(-0.15*t)';
ClusterFcn = (DipEvoFcn + 5).*B;
ModDepth = 1/ClusterFcn(1);
ClusterFcn = ClusterFcn/ClusterFcn(1);
ClusterFcn = ClusterFcn./sqrt(B);

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.005;
RegMatrix = regoperator(Dimension,3);
KB = dipolarkernel(t,r,ModDepth,sqrt(B));
Result = fitregmodel(ClusterFcn,KB,r,RegMatrix,'tv',RegParam,'Solver','fnnls');

err = any(abs(Result - Distribution)>1.5e-2);
maxerr = max(abs(Result - Distribution));

data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,Distribution,'k') 
    plot(r,Result,'r')
end

end