function [pass,maxerr] = test(opt)

% Check that apt() works with limited excitation bandwidth

t = linspace(0,3,200);
r = time2dist(t);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;
% Test apt using a 1GHz excitation bandwidth
aptK = aptkernel(t,'ExcitationBandwidth',1000);
[aptP,aptr] = apt(S,aptK,0.1);
P = dd_gauss(aptr,[3,0.5]);

% Pass: if the signal fits well enough
pass = all(abs(aptP - P) < 9e-1);
 
maxerr = max(abs(aptP - P));

if opt.Display
    figure(8),clf
    hold on
    plot(aptr,P,aptr,aptP)
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    legend('Truth','Fit')
    axis tight, grid on, box on
end


end