function [err,data] = test(opt,olddata)

%======================================================
% Control the pdsdata class construts properly
%======================================================

%Parameters
Offset = 1;
ModulationDepth = 0.35;
DecayRate = 0.0005;
Length = 200;
TimeStep = 16;
Cutoff = 70;

TimeAxis = linspace(0,Length*TimeStep,Length);
%Construct some dipolar evolution function from Fresnel integral
dipevo = 1-2*fresnels(TimeAxis*2*pi*1/(15^3));
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
classStructure = pdsdata();
%Introduce the data
classStructure.TimeAxis = TimeAxis;
classStructure.ExpData = clustersignal;
%And let the class prepare the time traces
classStructure = prepareFormFactor(classStructure,Cutoff);


clustersignal = clustersignal./clustersignal(1);
FormFactor = FormFactor./FormFactor(1);

%Check for errors
err(1) = classStructure.Length~=Length;
err(2) = classStructure.TimeStep~=TimeStep;
err(3) = any(abs(classStructure.FormFactor - olddata.FormFactor)>1e-10);
err(4) = any(abs(classStructure.DipEvoFcn - olddata.DipEvoFcn)>1e-10);
err(5) = any(abs(classStructure.ClusterSignal - olddata.ClusterSignal)>1e-10);
err(6) = abs(classStructure.ModDepth - olddata.ModDepth)/ModulationDepth>0.05;

err = any(err);

if ~err
  data.FormFactor = classStructure.FormFactor;
  data.DipEvoFcn = classStructure.DipEvoFcn;
  data.ClusterSignal = classStructure.ClusterSignal;
  data.ModDepth = classStructure.ModDepth;
else
  data = [];
end

end