function [pass,maxerr] = test(opt)

% Make sure the fractal-dimension homogeneous background is correct

t_us = 5; % us
c = 1e-5; % umol dm^-d
lam = 0.7;

d = 2;
B = bg_homfractal(t_us,[c d],lam);

Bref = 0.67388;

maxerr = abs((B-Bref)/Bref);

pass = maxerr<1e-5;

end
