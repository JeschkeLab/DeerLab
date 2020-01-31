function [err,data,maxerr] = test(opt,olddata)

N = 100;
r = linspace(1,5,N);

K = dipolarkernel(0.5,r);

err = any(size(K,1)>1);
maxerr = 0;
data = [];

end