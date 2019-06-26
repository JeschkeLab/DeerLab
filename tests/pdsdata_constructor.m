function [err,data] = test(opt,olddata)

%======================================================
% Control the pdsdata class construts properly
%======================================================

%Parameters
Offset = 2;
ModulationDepth = 0.35;
DecayRate = 0.0005;
Length = 200;
TimeStep = 16;

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
myClass = pdsdata('TimeAxis',TimeAxis,'ExpData',clustersignal);
%And let the class prepare the time traces
myClass = prepareFormFactor(myClass);


clustersignal = clustersignal./clustersignal(1);
FormFactor = FormFactor./FormFactor(1);

%Check for errors
err(1) = myClass.Length~=Length;
err(2) = myClass.TimeStep~=TimeStep;
err(3) = any(abs(myClass.FormFactor - olddata.FormFactor)>1e-10);
err(4) = any(abs(myClass.DipEvoFcn - olddata.DipEvoFcn)>1e-10);
err(5) = any(abs(myClass.ClusterFcn - olddata.ClusterSignal)>1e-10);
err(6) = abs(myClass.ModDepth - olddata.ModDepth)/ModulationDepth>0.05;

err = any(err);

if ~err
  data.FormFactor = myClass.FormFactor;
  data.DipEvoFcn = myClass.DipEvoFcn;
  data.ClusterSignal = myClass.ClusterFcn;
  data.ModDepth = myClass.ModDepth;
else
  data = [];
end

end