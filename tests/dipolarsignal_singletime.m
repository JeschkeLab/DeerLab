function [err,data,maxerr] = test(opt,olddata)

N = 100;
r = linspace(1,5,N);
P = rd_onegaussian(r,[3 0.5]);

V = dipolarsignal(0.5,r,P);

err = numel(V)>1;
maxerr = 0;
data = [];

end