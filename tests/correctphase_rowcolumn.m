function [err,data,maxerr] = test(opt,olddata)


V0 = 1:100;
inputPhase = pi/4;

V = V0.*exp(-1i*inputPhase);

V1 = correctphase(V);
V2 = correctphase(V.');

err = any(abs(V1 - V2)>1e-10);
maxerr = max(abs(V1 - V2));
data = [];

end