function [pass,maxerr] = test(opt)

% Check that fitparamodel's goodness of fit output works with global fitting

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

[~,Pfit,~,~,stats] = fitparamodel({V1,V2},@dd_gauss,r,{K1,K2});

Vfit1 = K1*Pfit;
Vfit2 = K2*Pfit;

dof = 2;
chi2red1 = 1/(numel(V1)-dof)*norm(V1 - Vfit1)^2/sigma^2;
chi2red2 = 1/(numel(V2)-dof)*norm(V2 - Vfit2)^2/sigma^2;


% Pass 1-2: the internal chi2red is computed correctly
pass(1) = abs(chi2red1 - stats{1}.chi2red) < 5e-2;
pass(2) = abs(chi2red2 - stats{2}.chi2red) < 5e-2;

maxerr = max(abs(chi2red1 - stats{1}.chi2red),abs(chi2red2 - stats{2}.chi2red));

end