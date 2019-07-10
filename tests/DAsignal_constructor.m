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
% dipevo = formfactor./bckg - ModulationDepth;
% dipevo = dipevo./dipevo(1);

%Cosntruct the class to be tested
myClass = DAsignal('TimeAxis',TimeAxis,'ExpData',clustersignal);
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

if opt.Display
figure(8),clf
subplot(131),hold on
plot(TimeAxis,FormFactor)
plot(TimeAxis,myClass.FormFactor)

subplot(132),hold on
plot(TimeAxis,dipevo)
plot(TimeAxis,myClass.DipEvoFcn)

subplot(133),hold on
plot(TimeAxis,clustersignal)
plot(TimeAxis,myClass.ClusterFcn)
end
    


end