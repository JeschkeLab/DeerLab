function [err,data,maxerr] = test(opt,olddata)

N = 100;
t = linspace(0,3,N);

K = dipolarkernel(t,3);

err = any(size(K,2)>1);
maxerr = 0;
data = [];

end