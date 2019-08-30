function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

clear apt 
%Parameters
Dimension = 300;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

DistDomainSmoothing = 0.2;
%Test apt using a 1GHz excitation bandwidth
aptKernel = aptkernel(TimeAxis,'ExcitationBandwidth',1000);
tic
[postDistribution] = apt(DipEvoFcn,aptKernel,DistDomainSmoothing);
pre = toc;
tic
[preDistribution] = apt(DipEvoFcn,aptKernel,DistDomainSmoothing);
post = toc;

error = abs(postDistribution - preDistribution);
err(1) = any(error>1e-20);
err(2) = post > pre;
data = [];
maxerr = max(error);

if opt.Display
figure(8),clf
hold on
plot(DistanceAxis,postDistribution)
plot(DistanceAxis,postDistribution)
end

end