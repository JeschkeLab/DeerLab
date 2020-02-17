function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,3,100);

K = dipolarkernel(t,3);

err = size(K,2)~=1;

maxerr = 0;
data = [];

end