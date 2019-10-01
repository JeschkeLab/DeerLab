function [err,data,maxerr] = test(opt,olddata)

N = 800;
dt = 0.008;
t = linspace(0,dt*(N-1),N);

tic 
preK = aptkernel(t);
pre = toc;
tic 
postK = aptkernel(t);
post = toc;

err(1) = post>pre/10;
err(2) = any(any(abs(preK.Base - postK.Base)>1e-15));
err = any(err);
data = [];
maxerr = max(max(abs(preK.Base - postK.Base)));


end