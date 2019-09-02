function [err,data] = test(opt,olddata)

N = 200;
dt = 0.008;
t = linspace(0,dt*(N-1),N);

%us
K1 = aptkernel(t);
%ns
t = t*1000;
K2 = aptkernel(t);

err = any(any(abs(K1.Base-K2.Base)>1e-10));
data = [];

end