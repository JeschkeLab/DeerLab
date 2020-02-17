function [err,data,maxerr] = test(opt,olddata)

N = 100;
t = linspace(0,3,N);
r = time2dist(t);

K1 = dipolarkernel(t,r,'Method','grid');
K2 = dipolarkernel(t,r,'Method','grid','Knots',1201);

delta = abs(K1 - K2);
err = any(any(delta>1e-5));
maxerr = max(max(delta));
data = [];

end