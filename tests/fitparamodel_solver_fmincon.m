function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.5];
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

InitialGuess = [2 0.1];
[FitDistribution,FitParam] = fitparamodel(DipEvoFcn,K,r,@rd_onegaussian,[],'Solver','fmincon');
err(1) = any(abs(FitDistribution - Distribution)>1e-5);
err(2) = any(abs(FitParam - InputParam)>1e-3);
err = any(err);

maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitDistribution,'r')
end

end