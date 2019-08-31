function [err,data,maxerr] = test(opt,oldata)

clear fitparamodel

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.5];
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

InitialGuess = [2 0.1];
tic
[preFitDistribution,FitParam] = fitparamodel(DipEvoFcn,K,r,@rd_onegaussian,InitialGuess);
pre = toc;
tic
[postFitDistribution,FitParam] = fitparamodel(DipEvoFcn,K,r,@rd_onegaussian,InitialGuess);
post = toc;


err(1) = any(abs(preFitDistribution - postFitDistribution)>1e-15);
err(2) = post > pre/4;
err = any(err);

maxerr = max(abs(preFitDistribution - postFitDistribution));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(r,preFitDistribution,'b')
   plot(r,postFitDistribution,'r')
end

end