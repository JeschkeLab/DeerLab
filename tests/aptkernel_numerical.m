function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check apt kernel is constructed properly
%======================================================

load('oldaptkernel512');

N = 512;
dt = 0.008;
t = linspace(0,dt*(N-1),N);
K = aptkernel(t);

FreqAxis = K.FreqAxis;
Base = K.Base;
Base = Base/mean(abs(diff(FreqAxis)));
NormFactor = K.NormalizationFactor;
NormFactor = NormFactor/(mean(abs(diff(FreqAxis))))^2;
t = K.t;
CrossTalk = K.Crosstalk;

err(1) = any(any(abs(Base - base)>1e-3));
err(2) = any(abs(NormFactor - tnorm)>1e-1);
err(3) = any(abs(FreqAxis - ny)>1e-3);
err(4) = any(abs(t - t)>1e-3);
err(5) = any(any(abs(CrossTalk - crosstalk)>1e-2));
err = any(err);
maxerr = max(max(abs(Base - base)));
data = [];

end