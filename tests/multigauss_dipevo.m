function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [4 0.2 4 1 3 0.4 0.4 0.4];
Distribution = rd_threegaussian(r,InputParam);
Distribution = Distribution/(1/sqrt(2*pi)*1/InputParam(2));
Distribution = Distribution/sum(Distribution)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

[FitDistribution,FitParam] = multigauss(DipEvoFcn,K,r,5,'aicc');
err = any(abs(FitDistribution - Distribution)>1e-5);


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