function [err,data,maxerr] = test(opt,olddata)

%==============================================================
% Get start of background fit ensure that integer is returned
%==============================================================
clear backgroundstart
%Parameters
DecayRate = 0.5;
Length = 200;
dt = 0.016;
t = linspace(0,Length*dt,Length);
%Construct some dipolar evolution function 
r = time2dist(t);
dipevo = dipolarkernel(t,r)*rd_onegaussian(r,[3,0.5]);
%Construct background
bckg = exp(-DecayRate*t).';
lambdaoriginal = 0.5;
%Account modulation depth for the offset=1
FormFactor = (1 - lambdaoriginal) + lambdaoriginal*dipevo;
S = FormFactor.*bckg;

%us
FitStartTime1 = backgroundstart(S,t,@td_exp);
%ns
t = t*1000;
FitStartTime2 = backgroundstart(S,t,@td_exp);

%Check for errors
err = abs(FitStartTime1 - FitStartTime2)>1e-10;
maxerr = max(abs(FitStartTime1 - FitStartTime2));
data = [];

if opt.Display
    figure(8),clf
    plot(t,bckg,t,Bfit)
end

end