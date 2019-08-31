function [err,data,maxerr] = test(opt,oldata)


Dimension = 300;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [2 0.4 3.5 0.3 5 0.3 0.3 0.3];
Distribution = rd_threerice(r,InputParam);
Distribution = Distribution/sum(Distribution)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

[FitDistribution] = fitparamodel(DipEvoFcn,K,r,@rd_threerice,0.75*InputParam,'solver','lsqnonlin');
err = any(abs(FitDistribution - Distribution)>1e-2);

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