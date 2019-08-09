function [err,data] = test(opt,olddata)

N = 800;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*(N-1),N);

clear aptkernel

tic 
preKernel = aptkernel(TimeAxis);
pre = toc;
tic 
postKernel = aptkernel(TimeAxis);
post = toc;

err(1) = post>pre/10;
err(2) = any(any(abs(preKernel.Base - postKernel.Base)>1e-15));
err = any(err);
data = [];
maxerr = max(max(abs(preKernel.Base - postKernel.Base)));


end