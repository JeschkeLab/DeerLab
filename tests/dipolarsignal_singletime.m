function [pass,maxerr] = test(opt)

% Check that one can generate signals with one single time-domain point

r = linspace(1,5,100);
P = rd_onegaussian(r,[3 0.5]);
V = dipolarsignal(0.5,r,P);

% Pass: a single point is returned
pass = numel(V)==1;

maxerr = NaN;

end