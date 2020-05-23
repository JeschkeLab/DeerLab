function [pass,maxerr] = test(opt)

% Check the absolute values returned by selregparam 

rng(1)
t = linspace(0,5,200);
r = linspace(2,5,100);
P = dd_gauss(r,[3,0.4]);
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.02,'rescale');

methods = {'lr','lc','cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
alphaopt1 = selregparam(S,K,r,'tikhonov',methods,'NonNegConstrained',false,'NoiseLevel',0.05,'search','grid');

% Values returned by the selregparam function (18/05/2020)
alpharef1 = [0.0005 0.0005 0.0631 0.0003 0.0003 0.0016 0.0003 0.0005 0.0003 0.0008 0.0251 0.0006 0.0003 0.0008];

% Pass: calculated values are equal to reference values (log10 scale)
threshold = 1e-2;
pass(1) = all(abs(alpharef1 - alphaopt1) <= threshold);

pass = all(pass);

maxerr = max(abs(log10(alpharef1) - log10(alphaopt1)));

end
