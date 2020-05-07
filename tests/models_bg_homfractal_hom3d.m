function [pass,maxerr] = test(opt)

% Make sure fractal homogeneous funtion with d=3 yields same result as dedicated
% 3D background function

t_us = 10;
c_uM = 100;
lam = 0.7;

B1 = bg_hom3d(t_us,c_uM,lam);
B2 = bg_homfractal(t_us,[c_uM 3],lam);

maxerr = abs(B1-B2)/B1;
pass = maxerr<1e-14;

end
