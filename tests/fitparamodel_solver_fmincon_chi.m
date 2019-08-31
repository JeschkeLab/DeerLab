function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [2 0.3];
Distribution = rd_onegaussian(r,[2,0.3]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

[FitDistribution] = fitparamodel(DipEvoFcn,K,r,@rd_onegaussian,0.6*InputParam,'Solver','fmincon','costmodel','chisquare');
err(1) = any(abs(FitDistribution - Distribution)>1e-7);
err = any(err);

maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
     figure(1),clf
   subplot(121)
   hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitDistribution,'r')
   subplot(122)
   hold on
   plot(r,Distribution,'b')
   plot(r,FitDistribution,'r')
   legend('truth','fit')
end

end