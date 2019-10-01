function [err,data,maxerr] = test(opt,olddata)

N = 200;
dt = 0.008;
t = linspace(0,dt*(N-1),N);

K1 = aptkernel(t,'ExcitationBandwidth',1e8);
K2 = aptkernel(t);

err = any(any(abs(K1.Base-K2.Base)>1e-10));
maxerr = max(max(abs(K1.Base-K2.Base)));
data = [];

end