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

%Check whether control statements on too smal start/end values work
backgroundstart(S,t,@td_exp,'RelSearchStart',0.00001,'RelSearchEnd',0.00005);
%Check other options
backgroundstart(S,t,@td_exp,'EndCutoffPos',10);

FitStartTime = backgroundstart(S,t,@td_exp,'RelSearchStart',0.05,'RelSearchEnd',0.8);
[Bfit] = fitbackground(S,t,@td_exp,FitStartTime);

%Check for errors
err = abs(bckg - Bfit)>1e-5;
maxerr = max(abs(bckg - Bfit));
data = [];

if opt.Display
    figure(8),clf
    plot(t,bckg,t,Bfit)
end

end