function [pass,maxerr] = test(opt)

N = 100;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
S = K*P;

[OptParam,Functionals,RegParams] = selregparam(S,K,r,'tikhonov',{'aic','gml','gcv'},'Search','grid');
%Accept testif all values are the same (should be as there is no noise)
err(1) = any(diff(OptParam) > 1e-2);
pass = all(err);
 
maxerr = max(diff(OptParam));

if opt.Display
    figure(8),clf
    subplot(131)
    plot(log(RegParams),log(Functionals{1}+ 1e5),'.')
    subplot(132)
    plot(log(RegParams),log(Functionals{2}),'.')
    subplot(133)
    plot(log(RegParams),log(Functionals{3}),'.')
end


end