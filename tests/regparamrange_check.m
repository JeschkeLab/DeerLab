function [pass,maxerr] = test(opt)

N = 200;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);

K = dipolarkernel(t,r);
L = regoperator(N,2);
lgRes = 0.1;
alpha = regparamrange(K,L,'Resolution',lgRes);

err = length(alpha)~=85;
maxerr = NaN;
 

end