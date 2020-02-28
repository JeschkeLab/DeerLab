function [pass,maxerr] = test(opt)

% Check indifference of diplarkernel() towards input dimensionality

t = linspace(-1,4,20);
r = linspace(1,6,50);
B = td_exp(t,0.3);
lam = 0.5;

K1 = dipolarkernel(t,r,lam,B,'Cache',false);
K2 = dipolarkernel(t.',r,lam,B,'Cache',false);
K3 = dipolarkernel(t,r.',lam,B,'Cache',false);
K4 = dipolarkernel(t,r,lam,B.','Cache',false);
K5 = dipolarkernel(t.',r.',lam,B,'Cache',false);
K6 = dipolarkernel(t,r.',lam,B.','Cache',false);
K7 = dipolarkernel(t.',r,lam,B.','Cache',false);
K8 = dipolarkernel(t.',r.',lam,B.','Cache',false);

% Pass: all kernels are equal
pass = isequal(K1,K2,K3,K4,K5,K5,K6,K7,K8);

maxerr = NaN;
 

end