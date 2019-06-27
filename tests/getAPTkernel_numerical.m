function [err,data] = test(opt,olddata)

%======================================================
% Check APT kernel is constructed properly
%======================================================

load(fullfile('comparison','oldAPTkernel512'));

N = 512;
TimeStep = 0.008;

[Base,NormFactor,FreqAxis,TimeAxis,CrossTalk] = getAPTkernel(N,TimeStep);

err(1) = any(any(abs(Base - base)>1e-3));
err(2) = any(abs(NormFactor - tnorm)>1e-1);
err(3) = any(abs(FreqAxis - ny)>1e-3);
err(4) = any(abs(TimeAxis - t)>1e-3);
err(5) = any(any(abs(CrossTalk - crosstalk)>1e-2));
err = any(err);
data = [];

end