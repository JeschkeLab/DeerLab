function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 200;
t = linspace(-0.6,2,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

aptK = aptkernel(t);

[aptP1] = apt(DipEvoFcn,aptK,0.05);
[aptP2] = apt(DipEvoFcn.',aptK,0.05);

error = abs(aptP1 - aptP2);
err = any(error>1e-10);
data = [];
maxerr = max(error);


end