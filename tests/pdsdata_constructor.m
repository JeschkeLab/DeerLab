function [err,data] = test(opt,olddata)

%======================================================
% Control the pdsdata class construts properly
%======================================================

%Parameters
Offset = 3;
ModulationDepth = 0.4;
DecayRate = 0.0005;
TimeStep = 16;
Length = 200;

TimeAxis = linspace(16,Length*TimeStep,Length);
%Construct some dipolar evolution function from Fresnel integral
dipevo = 1 - 2*fresnels(TimeAxis*2*pi*1/(15^3));
dipevo = dipevo(5:end);
TimeAxis = TimeAxis(5:end);
Length = length(TimeAxis);

%Construct background
bckg = exp(-DecayRate*TimeAxis);
%Account modulation depth for the offset=1
adaptedmodulationDepth = ModulationDepth*(1+Offset);
FormFactor = (1 - adaptedmodulationDepth) + adaptedmodulationDepth*dipevo;
FormFactor = FormFactor + Offset;
clustersignal = FormFactor.*bckg;
% dipevo = formfactor./bckg - ModulationDepth;
% dipevo = dipevo./dipevo(1);

%Cosntruct the class to be tested
myClass = pdsdata('TimeAxis',TimeAxis,'ExpData',clustersignal);
%And let the class prepare the time traces
myClass = prepare(myClass);

clustersignal = clustersignal./clustersignal(1);
FormFactor = FormFactor./FormFactor(1);

%Check for errors
err(1) = myClass.Length~=Length;
err(2) = myClass.TimeStep~=TimeStep;
err(3) = any(abs(myClass.FormFactor - FormFactor)>1e-2);
err(4) = any(abs(myClass.DipEvoFcn - dipevo)>1e-1);
err(5) = any(abs(myClass.ClusterFcn - clustersignal)>1e-2);

err = any(err);
data = [];

end