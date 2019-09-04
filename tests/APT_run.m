function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
P = P/sum(P)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

DistDomainSmoothing = 0.1;
%Test apt using a 1GHz excitation bandwidth
aptK = aptkernel(t,'ExcitationBandwidth',1000);
[aptP,aptr] = apt(DipEvoFcn,aptK,DistDomainSmoothing);

error = abs(aptP - P);
err = any(error>9e-1);
data = [];
maxerr = max(error);

if opt.Display
    figure(8),clf
    hold on
    plot(r,P)
    plot(aptr,aptP)
    
end

end