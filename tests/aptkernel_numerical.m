function [err,data] = test(opt,olddata)

%======================================================
% Check apt kernel is constructed properly
%======================================================

load(fullfile('comparison','oldaptkernel512'));

N = 512;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*(N-1),N);
Kernel = aptkernel(TimeAxis);


Base = Kernel.Base;
NormFactor = Kernel.NormalizationFactor;
FreqAxis = Kernel.FreqAxis;
TimeAxis = Kernel.TimeAxis;
CrossTalk = Kernel.Crosstalk;

err(1) = any(any(abs(Base - base)>1e-3));
err(2) = any(abs(NormFactor - tnorm)>1e-1);
err(3) = any(abs(FreqAxis - ny)>1e-3);
err(4) = any(abs(TimeAxis - t)>1e-3);
err(5) = any(any(abs(CrossTalk - crosstalk)>1e-2));
err = any(err);
data = [];

end