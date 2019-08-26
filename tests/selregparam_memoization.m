function [err,data,maxerr] = test(opt,data)

%Clear persistent variable
clear selregparam

Dimension = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = Kernel*Distribution;

RegParamSet = regparamrange(Kernel,RegMatrix);

tic
[preOptParam,Functionals,RegParamRange] = selregparam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,'tikhonov',{'aic','gcv','lr'});
precached = toc;

tic
[postOptParam,Functionals,RegParamRange] = selregparam(RegParamSet,DipEvoFcn,Kernel,RegMatrix,'tikhonov',{'aic','gcv','lr'});
postcached = toc;

err(1) = postcached>=precached/10;
error = abs(preOptParam - postOptParam);
err(2) = any(any(error>1e-18));

err = any(err);
data = [];
maxerr = max(max(error));

end