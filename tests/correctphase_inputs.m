function [pass,maxerr] = test(opt)

% Check that apt() works with defined input schemes

V0 = 1:100;
V0 = V0(:);
inputPhase = pi/4;
V = V0.*exp(-1i*inputPhase);

%Scheme 1
Vout1 = correctphase(V);
%Scheme 2
Vout2 = correctphase(V,inputPhase);
%Scheme 3
Vout3 = correctphase(V,inputPhase,false);

%Pass 1: input scheme 1 returns correct fit 
pass(1) = all(abs(Vout1 - V0) < 1e-10);
%Pass 2: input scheme 2 returns correct fit 
pass(2) = all(abs(Vout2 - V0) < 1e-10);
%Pass 3: input scheme 3 returns correct fit 
pass(3) = all(abs(Vout3 - V0) < 1e-10);

pass = all(pass);

maxerr = max(abs(Vout3 - V0));
 

end