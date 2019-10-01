function [err,data,maxerr] = test(opt,olddata)

N = 100;
t = linspace(0,3,N);
r = time2dist(t);

K1 = dipolarkernel(t,r,'Method','explicit');
K2 = dipolarkernel(t,r,'Method','explicit','Knots',1201);

err = any(any(abs(K1 - K2)>1e-5));
maxerr = max(max(abs(K1 - K2)));
data = [];

end