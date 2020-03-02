function [pass,maxerr] = test(opt)

% Test selregparam with Huber regularization

t = linspace(0,3,200);
r = linspace(2,6,100);
P = rd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

alphaopt = selregparam(S,K,r,'huber',{'aic','gcv'});

% Pass: both methods find the same solutions
pass = abs(diff(alphaopt)) < 1e-2;

maxerr = abs(diff(alphaopt));
 



end