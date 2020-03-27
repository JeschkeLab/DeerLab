function [pass,maxerr] = test(opt)

% Check that apt() works with negative times

t = linspace(-0.6,2,200);
r = time2dist(t);
P = dd_onegauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;
aptK = aptkernel(t);
[aptP,aptr] = apt(S,aptK,0.05);
P = dd_onegauss(aptr,[3,0.5]);

%Pass: if the signal fits well enough
pass = all(abs(aptP - P) < 9e-1);
 
maxerr = max(abs(aptP - P));

if opt.Display
    plot(aptr,P,aptr,aptP)
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    legend('Truth','Fit')
    axis tight, grid on, box on
end

end