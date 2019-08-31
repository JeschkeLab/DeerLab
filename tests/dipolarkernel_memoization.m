function [err,data,maxerr] = test(opt,data)

%Clear persistent variable
clear dipolarkernel

Dimension = 1000;
t = linspace(0,6,Dimension);
r = time2dist(t);

tic
preK = dipolarkernel(t,r);
precached = toc;

tic
postK = dipolarkernel(t,r);
postcached = toc;

err(1) = postcached>=precached/10;
error = abs(postK - preK);
err(2) = any(any(error>1e-18));
data = [];
err = any(err);
maxerr = max(max(error));

end