function [err,data,maxerr] = test(opt,olddata)

%==============================================================
% Get start of background fit ensure that integer is returned
%==============================================================
clear fitbackground
%Parameters
DecayRate = 0.25;
Length = 200;
TimeStep = 0.016;
TimeAxis = linspace(0,Length*TimeStep,Length);
%Construct some dipolar evolution function 
r = time2dist(TimeAxis);
dipevo = dipolarkernel(TimeAxis,r)*rd_onegaussian(r,[3,0.5]);
%Construct background
bckg = 1  - DecayRate*TimeAxis.';
lambdaoriginal = 0.5;
%Account modulation depth for the offset=1
FormFactor = (1 - lambdaoriginal) + lambdaoriginal*dipevo;
Signal = FormFactor.*bckg;

[FitStartTime] = backgroundstart(Signal,TimeAxis,@td_poly1);
[Bfit,lambda] = fitbackground(Signal,TimeAxis,@td_poly1,FitStartTime);
Bfit = Bfit*(1-lambda);
bckg = bckg*(1-lambdaoriginal);
%Check for errors
err = abs(bckg - Bfit)>1e-4;
maxerr = max(abs(bckg - Bfit));
data = [];

if opt.Display
    figure(8),clf
    plot(TimeAxis,bckg,TimeAxis,Bfit)
end

end