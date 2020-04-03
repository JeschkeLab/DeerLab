function [pass,maxerr] = test(opt)

% Check the absolute values returned by selregparam 

rng(1)
t = linspace(0,5,200);
r = linspace(2,5,100);
P = dd_gauss(r,[3,0.4]);
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.02);

methods = {'lr','lc','cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
alphaopt1 = selregparam(S,K,r,'tikhonov',methods,'NonNegConstrained',false,'NoiseLevel',0.05,'search','grid');
methods = {'cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
alphaopt2 = selregparam(S,K,r,'tikhonov',methods,'NonNegConstrained',false,'NoiseLevel',0.05,'search','golden');

% Values returned by the selregparam function at commit 8918a091ca350bb7fc04febff870e385f169a7f7   (10/03/2020)
alpharef1 = [0.5012 0.5012   79.4328 0.3162 0.3162 1.9953 0.3162 0.6310 0.3162 1.0000   25.1189 0.7943 0.3162 1.0000];
alpharef2 = [79.1674 0.3215 0.3355 1.7724 0.2975 0.5711 0.3445 0.9753   25.8199 0.7562 0.3475 0.9532];

% Pass 1: all values have the same absolute values
pass(1) = all(abs(alpharef1 - alphaopt1)./alpharef1 < 1e-3);
% Pass 2: all values have the same absolute values
pass(2) = all(abs(alpharef2 - alphaopt2)./alpharef2 < 1e-3);

pass = all(pass);

maxerr = max(abs(alpharef2 - alphaopt2)./alpharef2);

end