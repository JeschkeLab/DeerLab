function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [2.5 0.5 4 0.5 0.4];
Distribution = rd_twogaussian(r,InputParam);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;
B = exp(-0.15*t)';
ClusterFcn = (DipEvoFcn + 5).*B;
ModDepth = 1/ClusterFcn(1);
ClusterFcn = ClusterFcn/ClusterFcn(1);
ClusterFcn = ClusterFcn./sqrt(B);


KB = dipolarkernel(t,r,B,ModDepth,'KBType','sqrt');

InitialGuess = [2 0.1 5 0.1 0.1];
[FitDistribution] = fitparamodel(ClusterFcn,KB,r,@rd_twogaussian,InitialGuess);
err(1) = any(abs(FitDistribution - Distribution)>1e-5);
err = any(err);
maxerr = max(abs(FitDistribution - Distribution));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(t,ClusterFcn,'b')
   plot(t,KB*FitDistribution,'r')
   subplot(122)
   hold on
   plot(r,Distribution,'b')
   plot(r,FitDistribution,'r')
   
end

end