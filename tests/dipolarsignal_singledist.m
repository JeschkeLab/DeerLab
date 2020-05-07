function [pass,maxerr] = test(opt)

% Check that one can generate signals with one single distance-domain point

t = linspace(0,5,100);
r = 4;
V = dipolarsignal(t,r);

% Pass: same number of points as the time-axis are returned
pass = numel(V)==numel(t);

maxerr = NaN;

end