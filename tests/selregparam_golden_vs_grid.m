function [pass,maxerr] = test(opt)

% Check that the golden search algorithm is faster than the grid search

t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.05);

tic
[alpha1] = selregparam(S,K,r,'tikhonov','aic','Search','golden');
tictoc1 = toc;

tic
[alpha2] = selregparam(S,K,r,'tikhonov','aic','Search','grid');
tictoc2 = toc;

% Pass 1: both methods find similar results
pass(1) = abs(alpha1 - alpha2) < 1e-1;
% Pass 2: the golden search method is faster
pass(2) = tictoc1 < tictoc2;

pass = all(pass);
 
maxerr = max(abs(alpha1 - alpha2));


end