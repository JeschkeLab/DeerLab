function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 200;
t = linspace(-0.6,2,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

DistDomainSmoothing = 0.05;
%Test apt using a 1GHz excitation bandwidth
aptK = aptkernel(t,'ExcitationBandwidth',1000);
[aptP,aptr] = apt(DipEvoFcn,aptK,DistDomainSmoothing);

P = rd_onegaussian(aptr,[3,0.5]);
error = abs(aptP - P);
err = any(error>9e-1);
data = [];
maxerr = max(error);

if opt.Display
    figure(8),clf
    hold on
    plot(aptr,P)
    plot(aptr,aptP)
    
end

end