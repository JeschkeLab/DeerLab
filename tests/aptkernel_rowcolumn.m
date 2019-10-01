function [err,data,maxerr] = test(opt,olddata)

N = 200;
t = linspace(0,3,N);

%row
K2 = aptkernel(t);

%column
K1 = aptkernel(t.');

err = any(any(abs(K1.Base-K2.Base)>1e-10));
maxerr = max(max(abs(K1.Base-K2.Base)));
data = [];


end