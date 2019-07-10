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
myClass = DAsignal('TimeAxis',TimeAxis,'ExpData',clustersignal);
%And let the class prepare the time traces
myClass = prepare(myClass);

bckg = bckg./bckg(1);
bckgFitted = myClass.Background./myClass.Background(1);

%Check for errors
err = any(abs(bckgFitted - bckg)>1e-2);
data = [];

end