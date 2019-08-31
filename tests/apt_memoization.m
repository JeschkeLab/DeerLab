function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

clear apt 
%Parameters
Dimension = 300;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);
Distribution = Distribution/sum(Distribution);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

DistDomainSmoothing = 0.2;
%Test apt using a 1GHz excitation bandwidth
aptK = aptkernel(t,'ExcitationBandwidth',1000);
tic
[postDistribution] = apt(DipEvoFcn,aptK,DistDomainSmoothing);
pre = toc;
tic
[preDistribution] = apt(DipEvoFcn,aptK,DistDomainSmoothing);
post = toc;

error = abs(postDistribution - preDistribution);
err(1) = any(error>1e-20);
err(2) = post > pre;
data = [];
maxerr = max(error);

if opt.Display
figure(8),clf
hold on
plot(r,postDistribution)
plot(r,postDistribution)
end

end