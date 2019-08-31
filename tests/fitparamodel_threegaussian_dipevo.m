function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.3 5 0.3 5.5 0.5 0.3 0.2];
Distribution = rd_threegaussian(r,InputParam);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

[FitDistribution] = fitparamodel(DipEvoFcn,K,r,@rd_threegaussian,[],'solver','lsqnonlin');
err(1) = any(abs(FitDistribution - Distribution)>1e-1);
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
   
end

end