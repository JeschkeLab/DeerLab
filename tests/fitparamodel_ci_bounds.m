function [pass,maxerr] = test(opt)

% Test that confidence levels are constrained correctly
t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

lam = 0.2;
K = dipolarkernel(t,r,lam);
V = K*P + whitegaussnoise(t,0.1);


Vmodel  = @(t,par) 1 - par(1) + par(1)*dipolarkernel(t,r)*dd_gauss(r,par(2:3));
par0 = [0.5 4 0.5];
upper = [0.3 20 0.6];
lower = [0.1 1 0.1];
[~,~,parCI] = fitparamodel(V,Vmodel,t,par0,'lower',lower,'upper',upper);


% Pass 1-2: confidence intervals are within the bounds
pass(1) = all(all(parCI <= upper.'));
pass(2) = all(all(parCI >= lower.'));

maxerr = NaN;

end