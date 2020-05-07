function [pass,maxerr] = test(opt)

% Check indifference of fitregmodel() towards input dimensionality

t = linspace(0,5,50);
r = linspace(1,6,50);
P = dd_gauss(r,[4 0.3]);
K = dipolarkernel(t,r);
S = K*P;

Pfit1 = fitregmodel(S,K,r);
Pfit2 = fitregmodel(S.',K,r);
Pfit3 = fitregmodel(S,K,r.');
Pfit4 = fitregmodel(S.',K,r.');

% Pass 1: all distributions are equal
pass(1) = isequal(Pfit1,Pfit2,Pfit3,Pfit4);
% Pass 2: all distributions are column vectors
pass(2) = iscolumn(Pfit1) & iscolumn(Pfit2) & iscolumn(Pfit3) & iscolumn(Pfit4);

pass = all(pass);

maxerr = NaN;
 


end