function [pass,maxerr] = test(opt)

% Chech that number of knots in dipolarkernel() can be modified

t = linspace(0,3,30);
r = time2dist(t);

K1 = dipolarkernel(t,r,'Method','grid');
K2 = dipolarkernel(t,r,'Method','grid','nKnots',1201);

% Pass: both kernels are equal 
delta = abs(K1 - K2);
pass = all(delta(:) < 1e-4);

maxerr = max(delta(:));
 

end