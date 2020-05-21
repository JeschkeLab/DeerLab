function [pass,maxerr] = test(opt)

% Check that fitregmodel's goodness of fit output works with global fitting

rng(1)

t1 = linspace(0,8,500);
t2 = linspace(0,3,200);
r = linspace(2,5,300);
P = dd_gauss(r,[4.5 0.5]);
K1 = dipolarkernel(t1,r);
K2 = dipolarkernel(t2,r);

sigma = 0.01;
V1 = K1*P + whitegaussnoise(t1,sigma);
V2 = K2*P + whitegaussnoise(t2,sigma);

[~,~,~,stats] = fitregmodel({V1,V2},{K1,K2},r);

% Pass 1-2: the internal chi2red is computed correctly
pass(1) = abs(stats{1}.chi2red -1) < 1e-1;
pass(2) = abs(stats{2}.chi2red -1) < 1e-1;

maxerr = NaN;

end