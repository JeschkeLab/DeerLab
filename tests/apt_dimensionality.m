function [pass,maxerr] = test(opt)

% Check indifference of apt() towards input dimensionality

t = linspace(-1,4,100);
S = dipolarsignal(t,3);
K = aptkernel(t);
Pfit1 = apt(S,K);
Pfit2 = apt(S.',K);

% Pass 1: outputs are unchanged
pass(1) = isequal(Pfit1,Pfit2);
% Pass 2: both outputs are columns
pass(2) = iscolumn(Pfit1) & iscolumn(Pfit2);
pass = all(pass);

maxerr = NaN;

end