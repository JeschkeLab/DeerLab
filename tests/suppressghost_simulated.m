function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Test supressghost.m
%=======================================


Dimension = 100;
t = linspace(0.006,3,Dimension);
r = time2dist(t);

P = rd_onegaussian(r,[3,0.5]);
DistrAB = P/sum(P)';
P = rd_onegaussian(r,[4,0.5]);
DistrAC = P/sum(P)';
P = rd_onegaussian(r,[5,0.5]);
DistrBC = P/sum(P)';
P = DistrAB + DistrBC + DistrAC;
DistrAB = DistrAB/sum(P)';
DistrBC = DistrBC/sum(P)';
DistrAC = DistrAC/sum(P)';

FreqAxis = 52.04./(r.^3)';

wddt = 2*pi*FreqAxis.*t;
kappa = sqrt(6*wddt/pi);
%Compute Fresnel integrals of 0th order

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private
C = fresnelC(kappa);
S = fresnelS(kappa);
cd(currentpath)
omegaAB = sum(sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S).*DistrAB);
omegaBC = sum(sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S).*DistrBC);
omegaAC = sum(sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S).*DistrAC);

lambda = 0.6;

P = 1/3*(omegaAB + omegaBC + omegaAC);
T = 1/6*(omegaAB.*omegaAC + omegaAB.*omegaBC + omegaBC.*omegaAC);


FormFactor3 = 1-2*lambda + lambda^2 +lambda*(1-lambda)*P + lambda^2*T;
FormFactor2 = 1-lambda + lambda*P;

correctedFormFactor3 = FormFactor3.^(1/2);


Dip3 = (correctedFormFactor3-(1-lambda))/lambda;
Dip2 = (FormFactor2-(1-lambda))/lambda;

Dip2 = Dip2/Dip2(1);
Dip3 = Dip3/Dip3(1);

signal = suppressghost(FormFactor3,3);
signal = signal/signal(1);
signal = (signal - lambda)/(1-lambda) ;
signal= signal/signal(1);
K = dipolarkernel(t,r);
% K = K/mean(diff(r));
RegMatrix = regoperator(Dimension,2);
RegParam = 4;
Dist2 = fitregmodel(Dip2,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
Dist3 = fitregmodel(Dip3,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
DistrTest = fitregmodel(signal,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');

err(1) = any(abs(DistrTest - Dist3)>3e-1);
err(2) = any(abs(DistrTest - Dist2)>3e-1);

maxerr = max(abs(DistrTest - Dist2));
err = any(err);
data = [];

if opt.Display
    figure(8),clf
    hold on
    plot(r,Dist2,'k')
    plot(r,DistrTest,'r')
    
end



end

