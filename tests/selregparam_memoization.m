function [err,data,maxerr] = test(opt,data)

%Clear persistent variable
clear selregparam

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);
Distribution = Distribution/sum(Distribution);

K = dipolarkernel(t,r);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = K*Distribution;

RegParamSet = regparamrange(K,RegMatrix);

tic
[preOptParam,Functionals,RegParamRange] = selregparam(RegParamSet,DipEvoFcn,K,RegMatrix,'tikhonov',{'aic','gcv','lr'});
precached = toc;

tic
[postOptParam,Functionals,RegParamRange] = selregparam(RegParamSet,DipEvoFcn,K,RegMatrix,'tikhonov',{'aic','gcv','lr'});
postcached = toc;

err(1) = postcached>=precached/10;
error = abs(preOptParam - postOptParam);
err(2) = any(any(error>1e-18));

err = any(err);
data = [];
maxerr = max(max(error));

end