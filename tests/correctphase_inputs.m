function [err,data,maxerr] = test(opt,olddata)

V0 = 1:100;
inputPhase = pi/4;

V = V0.*exp(-1i*inputPhase);

Vout1 = correctphase(V);
Vout2 = correctphase(V,inputPhase);
Vout3 = correctphase(V,inputPhase,false);
Vout4 = correctphase(V,[],true);

err(1) = any(abs(Vout1 - V0)>1e-5);
err(2) = any(abs(Vout2 - V0)>1e-5);
err(3) = any(abs(Vout3 - V0)>1e-5);
err(4) = any(abs(Vout4 - V0)>1e-2);

err = any(err);
maxerr = max(abs(Vout1 - V0));
data = [];

end