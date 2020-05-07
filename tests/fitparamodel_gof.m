function [pass,maxerr] = test(opt)

% Check that fitparamodel's goodness of fit output works

rng(1)

t = linspace(0,8,500);
r = linspace(2,5,300);
P = dd_gauss(r,[4.5 0.5]);
K = dipolarkernel(t,r);

sigma = 0.01;
V = K*P + whitegaussnoise(t,sigma);

[~,Pfit,~,~,stats] = fitparamodel(V,@dd_gauss,r,K);
Vfit = K*Pfit;

dof = 2;
chi2red = 1/(numel(V)-dof)*norm(V - Vfit)^2/sigma^2;

% Pass: the internal chi2red is computed correctly
pass = abs(chi2red - stats.chi2red) < 3e-2;

maxerr = max(abs(chi2red - stats.chi2red));

end