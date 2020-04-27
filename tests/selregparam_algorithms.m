function [pass,maxerr] = test(opt)

% Check that the all search algorithms find the same minimum

rng(553)
t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
V = dipolarsignal(t,r,P,'noiselevel',0.05);

RegMethod = 'Tikhonov';
SelMethod = 'AIC';

[alpha1] = selregparam(V,K,RegMethod,SelMethod,'Search','fminbnd');
[alpha2] = selregparam(V,K,RegMethod,SelMethod,'Search','golden');
[alpha3,~,alphas] = selregparam(V,K,RegMethod,SelMethod,'Search','grid');

lga = log10([alpha1 alpha2 alpha3]);

% Pass if all methods find similar results
pass(1) = abs(lga(2)-lga(1)) < 0.02;

dlga = abs(mean(diff(log10(alphas))));
pass(2) = abs(lga(3)-lga(1)) < dlga/2;
 
maxerr = max(abs(alpha2 - alpha3));

end