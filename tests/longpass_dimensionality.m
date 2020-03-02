function [pass,maxerr] = test(opt)

% Check indifference of longpass() towards input dimensionality

t = linspace(-1,4,100);
S = dipolarsignal(t,3);

S1 = longpass(t,S);
S2 = longpass(t.',S);
S3 = longpass(t,S.');
S4 = longpass(t.',S.');

% Pass 1: all signals are equal
pass(1) = isequal(S1,S2,S3,S4);
% Pass 2: all signals are column vectors
pass(2) = iscolumn(S1) & iscolumn(S2) & iscolumn(S3) & iscolumn(S4);

pass = all(pass);

maxerr = NaN;
 

end