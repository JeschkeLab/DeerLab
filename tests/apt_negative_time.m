function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 200;
dt = 0.008;
t = linspace(-0.6,2,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);
Distribution = Distribution/sum(Distribution)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

DistDomainSmoothing = 0.05;
%Test apt using a 1GHz excitation bandwidth
aptK = aptkernel(t,'ExcitationBandwidth',1000);
[aptDistribution,aptr] = apt(DipEvoFcn,aptK,DistDomainSmoothing);

error = abs(aptDistribution - Distribution);
err = any(error>9e-1);
data = [];
maxerr = max(error);

if opt.Display
    figure(8),clf
    subplot(121)
    hold on
    plot(t,DipEvoFcn)
    plot(t,dipolarkernel(t,aptr)*aptDistribution)
    subplot(122)
    hold on
    plot(r,Distribution)
    plot(aptr,aptDistribution)
    
end

end