function [pass,maxerr] = test(opt)

% Test that multiple confidence levels can be requested via options

t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.05);

[~,~,parfit,parCI] = selectmodel({@dd_gauss,@dd_gauss2},S,r,K,'aic','confidencelevel',[0.5 0.95]);

parci1 = parCI{1}{1};
parci2 = parCI{1}{2};

% Pass 1-2: confidence intervals behave as expected
pass = all(all(abs(parfit{1} - parci1.') < abs(parfit{1} - parci2.')));

maxerr = NaN;


end