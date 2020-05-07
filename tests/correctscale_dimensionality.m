function [pass,maxerr] = test(opt)

% Check indifference of correctscale() towards input dimensionality

t = linspace(-1,4,100);
S = dipolarsignal(t,3,'scale',1e4);
t = t + abs(min(t));

V1 = correctscale(S,t);
V2 = correctscale(S.',t);
V3 = correctscale(S,t.');
V4 = correctscale(S.',t.');

% Pass 1: all signals are equal
pass(1) = isequal(V1,V2,V3,V4);
% Pass 2: all signals are columns
pass(2) = iscolumn(V1) & iscolumn(V2) & iscolumn(V3) & iscolumn(V4);

pass = all(pass);

maxerr = NaN;
 

end