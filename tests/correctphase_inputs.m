function [pass,maxerr] = test(opt)

% Check that various argument combinations work

V0 = 1:100;
V0 = V0(:);
ph = pi/4;
V = V0.*exp(1i*ph);

Vout1 = correctphase(V);
Vout2 = correctphase(V,ph);
Vout3 = correctphase(V,ph,false);

thr = 1e-10;
pass(1) = all(abs(Vout1 - V0) < thr);
pass(2) = all(abs(Vout2 - V0) < thr);
pass(3) = all(abs(Vout3 - V0) < thr);

pass = all(pass);

maxerr = max(abs(Vout3 - V0));
 

end