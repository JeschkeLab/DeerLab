function [pass,maxerr] = test(opt)

% Test selregparam with unconstrained Huber regularization

t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

[alphaopt,fcns,alphas] = selregparam(S,K,r,'huber',{'aic','gcv'},'NonNegConstrained',false);

% Pass: both methods find the same solutions
pass = abs(diff(alphaopt)) < 1e-2;

maxerr = abs(diff(alphaopt));
 

if opt.Display
   plot(log10(alphas),fcns{1},log10(alphas),fcns{2})
   xlabel('log_{10}(\alpha)')
   ylabel('Functional')
   legend('AIC','BIC')
   grid on, axis tight, box on
end


end