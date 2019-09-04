function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Dimension = 200;
dt = 0.008;
t = linspace(0,Dimension*dt,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
P = P/sum(P);
ExcitationBandwidth = 0.05; %MHz

K = dipolarkernel(t,r,'ExcitationBandwidth',ExcitationBandwidth);

Trace  = K*P;

err = any(abs(Trace - 0*Trace)>1e-5);
maxerr = max(abs(Trace -  0*Trace));
data = [];
end