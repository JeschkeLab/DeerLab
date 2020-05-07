function [pass,maxerr] = test(opt)

% Check the absolute values returned by selregparam 

rng(1)
t = linspace(0,5,200);
r = linspace(2,5,100);
P = dd_gauss(r,[3,0.4]);
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.02,'rescale');

methods = {'lr','lc','cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
alphaopt1 = selregparam(S,K,'tikhonov',methods,'NonNegConstrained',false,'NoiseLevel',0.05,'search','grid');

% Values returned by the selregparam function at commit cec7a2b5af55a10739ec7993dd1aee4c25920f6c  (13/04/2020)
alpharef1 = [0.5012 0.5012   79.4328 0.3162 0.3162 1.5849 0.3162 0.6310 0.3162 1.0000 25.1189 0.7943 0.3162 1.0000];

% Pass: calculated values are equal to reference values (log10 scale)
threshold = 1e-2;
pass(1) = all(abs(log10(alpharef1) - log10(alphaopt1)) <= threshold);

pass = all(pass);

maxerr = max(abs(log10(alpharef1) - log10(alphaopt1)));

end
