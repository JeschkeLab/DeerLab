function [err,data,maxerr] = test(opt,oldata)

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.5];
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;
B = exp(-0.15*t)';
ClusterFcn = (DipEvoFcn + 5).*B;
ModDepth = 1/ClusterFcn(1);
ClusterFcn = ClusterFcn/ClusterFcn(1);
ClusterFcn = ClusterFcn./sqrt(B);

KB = dipolarkernel(t,r,ModDepth,sqrt(B));

InitialGuess = [2 0.1];
[FitDistribution,FitParam] = fitparamodel(ClusterFcn,KB,r,@rd_onegaussian,InitialGuess);
err(1) = any(abs(FitDistribution - Distribution)>1e-5);
err(2) = any(abs(FitParam - InputParam)>1e-3);
err = any(err);

maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,ClusterFcn,'b')
   plot(t,KB*FitDistribution,'r')
end

end