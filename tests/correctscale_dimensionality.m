function [pass,maxerr] = test(opt)

% Check indifference of correctscale() towards input dimensionality

t = linspace(-1,4,100);
V = dipolarsignal(t,3,'scale',1e4);

V1 = correctscale(V,t);
V2 = correctscale(V.',t);
V3 = correctscale(V,t.');
V4 = correctscale(V.',t.');

% Pass 1: all signals are equal
pass(1) = isequal(V1,V2,V3,V4);
% Pass 2: all signals are columns
pass(2) = iscolumn(V1) & iscolumn(V2) & iscolumn(V3) & iscolumn(V4);

pass = all(pass);

maxerr = NaN;
 

end