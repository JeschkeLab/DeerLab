function [err,data,maxerr] = test(opt,oldata)


Ntime = 100;
Ndist = 200;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
r = linspace(rmin,rmax,Ndist);
InputParam = [3,0.5];
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

InitialGuess = [2 0.1];
[FitDistribution,FitParam] = fitparamodel(DipEvoFcn,K,r,@rd_onegaussian,InitialGuess);
err(1) = any(abs(FitDistribution - Distribution)>1e-5);
err(2) = any(abs(FitParam - InputParam)>1e-3);
err(3)  = length(FitDistribution) < length(DipEvoFcn);
err = any(err);

maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitDistribution,'r')
end

end