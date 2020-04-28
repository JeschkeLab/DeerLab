function [pass,maxerr] = test(opt)

% Make sure the 3D homogeneous background is correct

t_us = 5; % us
c_uM = 100; % uM
lam = 0.7;

B = bg_hom3d(t_us,c_uM,lam);

kapprox = 1e-3; % approximate decay constant, uM^-1 us^-1
Bapprox = exp(-kapprox*c_uM*lam*t_us);

maxerr = abs((B-Bapprox)/Bapprox);

pass = maxerr<1e-3;

end
