function [err,data] = test(opt,olddata)

%======================================================
% Control the pdsdata class construts properly
%======================================================

%Parameters
Offset = 1;
ModulationDepth = 0.35;
DecayRate = 0.0001;
tmax = 6032;
TimeStep = 8;
Cutoff = 200;

TimeAxis = 0:TimeStep:tmax;
%Construct some dipolar evolution function from Fresnel integral
dipevo = 1-2*fresnels(TimeAxis*2*pi*1/(19^3));
%Construct background
bckg = exp(-DecayRate*TimeAxis);
%Account modulation depth for the offset=1
adaptedmodulationDepth = ModulationDepth*(1+Offset);
formfactor = (1 - adaptedmodulationDepth) + adaptedmodulationDepth*dipevo;
formfactor = formfactor + Offset;
clustersignal = formfactor.*bckg;
dipevo = formfactor./bckg - adaptedmodulationDepth;
dipevo = dipevo./dipevo(1);

%Cosntruct the class to be tested
classStructure = pdsdata();
%Introduce the data
classStructure.TimeAxis = TimeAxis;
classStructure.ExpData = clustersignal;
%And let the class prepare the time traces 
classStructure = prepareFormFactor(classStructure,Cutoff);


clustersignal = clustersignal./clustersignal(1);

%Check for errors
err = classStructure.Length~=755;
err = classStructure.TimeStep~=TimeStep;
err = classStructure.TimeStep~=8;
err = any(abs(classStructure.FormFactor - formfactor)>1e-10);
err = any(abs(classStructure.DipEvoFcn - dipevo)>1e-10);
err = any(abs(classStructure.ClusterSignal - clustersignal)>1e-10);
err = abs(classStructure.ModDepth - ModulationDepth)/ModulationDepth>0.05;

data = [];

end