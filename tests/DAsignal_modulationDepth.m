function [err,data] = test(opt,olddata)

%======================================================
% Control the pdsdata class construts properly
%======================================================

%Parameters
Offset = 3;
ModulationDepth = 0.4;
DecayRate = 0.0005;
Length = 200;
TimeStep = 16;
TimeAxis = linspace(16,Length*TimeStep,Length);
%Construct some dipolar evolution function from Fresnel integral
rmin = (4*TimeStep/1000*52.04/0.85)^(1/3);
rmax = 6*(Length*TimeStep/1000/2)^(1/3);
TimeAxis = linspace(0,TimeStep*Length,Length);
Length = length(TimeAxis);

DistanceAxis = linspace(rmin,rmax,Length);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);
Kernel = getKernel(Length,TimeStep,rmin,rmax);

dipevo = Kernel*Distribution;
dipevo = dipevo';
%Construct background
bckg = exp(-DecayRate*TimeAxis);
%Account modulation depth for the offset=1
adaptedmodulationDepth = ModulationDepth*(1+Offset);
FormFactor = (1 - adaptedmodulationDepth) + adaptedmodulationDepth*dipevo;
FormFactor = FormFactor + Offset;
clustersignal = FormFactor.*bckg;


%Cosntruct the class to be tested
myClass = DAsignal('TimeAxis',TimeAxis,'ExpData',clustersignal);
%And let the class prepare the time traces
myClass = prepare(myClass);

%Check for errors
err = abs(myClass.ModDepth - ModulationDepth)>1e-2;
data = [];
end