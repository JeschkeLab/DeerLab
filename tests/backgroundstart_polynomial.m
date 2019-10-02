function [err,data,maxerr] = test(opt,olddata)

%==============================================================
% Get start of background fit ensure that integer is returned
%==============================================================

%Parameters
k = 0.25;
N = 200;
dt = 0.016;
t = linspace(0,N*dt,N);
%Construct some dipolar evolution function 
r = time2dist(t);
S = dipolarkernel(t,r)*rd_onegaussian(r,[3,0.5]);
%Construct background
bckg = 1  - k*t.';
lam0 = 0.5;
%Account modulation depth for the offset=1
S = (1 - lam0) + lam0*S;
S = S.*bckg;

[FitStartTime] = backgroundstart(S,t,@td_poly1);
[Bfit,lam] = fitbackground(S,t,@td_poly1,FitStartTime);
Bfit = Bfit*(1-lam);
bckg = bckg*(1-lam0);
%Check for errors
err = abs(bckg - Bfit)>1e-4;
maxerr = max(abs(bckg - Bfit));
data = [];

if opt.Display
    figure(8),clf
    plot(t,bckg,t,Bfit)
end

end