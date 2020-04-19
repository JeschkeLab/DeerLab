function [pass,maxerr] = test(opt)

% Check that selregparam returns the proper vectors of evaluated alphas

t = linspace(0,5,100);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.5]);
K = dipolarkernel(t,r);
V = K*P + whitegaussnoise(t,0.005);

[alphaopt1,aic1,alphas1] = selregparam(V,K,r,'tikh','aic','Search','golden');
[alphaopt2,aic2,alphas2] = selregparam(V,K,r,'tikh','aic','Search','grid');

% Pass 1-2: vectors have the right size
pass(1) = numel(aic1)==numel(alphas1);
pass(2) = numel(aic2)==numel(alphas2);

% Pass 1-2: vectors have the right size
pass(3) = iscolumn(alphas1);
pass(4) = iscolumn(alphas2);

% Pass 3-4: the optimal alpha is contained in the evaluated alphas
pass(5) = any(alphaopt1 == alphas1);
pass(6) = any(alphaopt2 == alphas2);

pass = all(pass);

maxerr = NaN;

end
