function [pass,maxerr] = test(opt)

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);

K = dipolarkernel(t,r);
L = regoperator(Dimension,2);

alpha = regparamrange(K,sparse(L));

err = (length(alpha)~=85);
maxerr = length(alpha)-85;
 

end