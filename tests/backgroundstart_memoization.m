function [err,data,maxerr] = test(opt,olddata)

%==============================================================
% Get start of background fit ensure that integer is returned
%==============================================================
clear backgroundstart
%Parameters
Offset = 2;
ModulationDepth = 0.35;
DecayRate = 0.0005;
Length = 200;
TimeStep = 16;
TimeAxis = linspace(TimeStep,Length*TimeStep,Length) - TimeStep;
%Construct some dipolar evolution function from Fresnel integral

r = time2dist(TimeAxis);
dipevo = dipolarkernel(TimeAxis,r)*rd_onegaussian(r,[3,0.5]);
%Construct background
bckg = exp(-DecayRate*TimeAxis);
%Account modulation depth for the offset=1
adaptedmodulationDepth = ModulationDepth*(1+Offset);
FormFactor = (1 - adaptedmodulationDepth) + adaptedmodulationDepth*dipevo;
FormFactor = FormFactor + Offset;
Signal = FormFactor.*bckg;

tic
preFitStartTime = backgroundstart(Signal,TimeAxis,@td_exp);
pre = toc;
tic
postFitStartTime = backgroundstart(Signal,TimeAxis,@td_exp);
post = toc;


%Check for errors
err(1) = any(abs(preFitStartTime - postFitStartTime)>1e-15);
err(2) = post > pre/4;

maxerr = max(abs(preFitStartTime - postFitStartTime));
err = any(err);
data = [];

end