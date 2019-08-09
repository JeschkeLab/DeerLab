function [err,data] = test(opt,olddata)

%==============================================================
% Get start of background fit ensure that integer is returned
%==============================================================

%Parameters
Offset = 2;
ModulationDepth = 0.35;
DecayRate = 0.0005;
Length = 200;
TimeStep = 16;
TimeAxis = linspace(TimeStep,Length*TimeStep,Length) - TimeStep;
%Construct some dipolar evolution function from Fresnel integral
dipevo = 1-2*fresnelS(TimeAxis*2*pi*1/(15^3));
%Construct background
bckg = exp(-DecayRate*TimeAxis);
%Account modulation depth for the offset=1
adaptedmodulationDepth = ModulationDepth*(1+Offset);
FormFactor = (1 - adaptedmodulationDepth) + adaptedmodulationDepth*dipevo;
FormFactor = FormFactor + Offset;
Signal = FormFactor.*bckg;


FitStartTime = backgroundstart(Signal,TimeAxis);

%Check for errors
err = abs(FitStartTime - 320)>0;
data = [];

end