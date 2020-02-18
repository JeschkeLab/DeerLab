function [err,data,maxerr] = test(opt,olddata)

N = 30;
t = linspace(0,3,N);
r = time2dist(t);

K1 = dipolarkernel(t,r,'Method','grid');
K2 = dipolarkernel(t,r,'Method','grid','Knots',1201);

delta = abs(K1 - K2);
err = any(delta(:)>1e-5);
maxerr = max(delta(:));
data = [];

end