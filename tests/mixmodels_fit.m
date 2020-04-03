function [pass,maxerr] = test(opt)

% Check that the mixed models can be used for fitting

t = linspace(0,3,200);
r = linspace(2,6,100);
parIn1 = [3 0.5];
parIn2 = [4 0.5];
mixedModel = mixmodels({@dd_gauss,@dd_gauss});
parInMix = [0.3 parIn1 parIn2];
Pmix = mixedModel(r,parInMix);

K = dipolarkernel(t,r);
S = K*Pmix;
par0 = 0.9*parInMix;
[~,Pfit] = fitparamodel(S,mixedModel,r,K,par0);

% Pass: the distribution is well fitted with the mixed model
pass = all(abs(Pmix - Pfit) < 1e-5);

maxerr = max(abs(Pmix - Pfit));
 
if opt.Display
   plot(r,Pmix,'k',r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end