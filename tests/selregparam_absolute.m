function [err,data,maxerr] = test(opt,olddata)

% Test the absolute values returned by selregparam 
% Fails if the values are different than on the 22/02/2020

t = linspace(0,5,200);
r = linspace(2,5,100);
P = rd_onegaussian(r,[3,0.4]);
P = P/sum(P);

K = dipolarkernel(t,r);
rng(1)
S = K*P + whitegaussnoise(t,0.02);

methods = {'lr','lc','cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
alphaopt1 = selregparam(S,K,r,'tikhonov',methods,'NonNegConstrained',false,'NoiseLevel',0.05,'search','grid');
methods = {'cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
alphaopt2 = selregparam(S,K,r,'tikhonov',methods,'NonNegConstrained',false,'NoiseLevel',0.05,'search','golden');

%Values returned by the selregparam function at commit 6a01f84addb76e0d4a4832ba206ede80b934aad2 
logalpha1 = [-5.7565 -6.6775 -11.5129 2.9934 2.9934 4.3749 2.9934 4.3749 2.9934 4.3749 4.3749 4.3749 3.6841 4.3749];
logalpha2 = [-11.4187 2.9825 3.0426 4.3148 2.9825 4.3148 3.0426 4.3148 4.3148 4.3148 3.6116 4.3148];

%Accept testif all values are the same (should be as there is no noise)
err(1) = any(abs(logalpha1 - log(alphaopt1))>1e-2);
err(2) = any(abs(logalpha2 - log(alphaopt2))>1e-2);
err = any(err);

maxerr = max(abs(logalpha2 - log(alphaopt2)));
data = [];



end