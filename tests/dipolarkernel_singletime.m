function [err,data,maxerr] = test(opt,olddata)

r = linspace(1,5,100);

K = dipolarkernel(0.5,r);

err = size(K,1)~=1;

maxerr = 0;
data = [];

end