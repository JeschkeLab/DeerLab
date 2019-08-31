function [err,data,maxerr] = test(opt,data)

%Clear persistent variable
clear regularize

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = K*Distribution;

RegParam = 100;

tic
preDist = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
precached = toc;

tic
postDist = fitregmodel(DipEvoFcn,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
postcached = toc;

err(1) = postcached>=precached/4;
error = abs(preDist - postDist);
err(2) = any(any(error>1e-18));
err = any(err);
data = [];
maxerr = max(max(error));

end