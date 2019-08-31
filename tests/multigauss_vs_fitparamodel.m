function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [2 0.5 3 0.5 0.4];
Distribution = rd_twogaussian(r,InputParam);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

InitialGuess = [2 0.1 5 0.1 0.5];
[FitDistribution,FitParam] = fitparamodel(DipEvoFcn,K,r,@rd_twogaussian,InitialGuess);
[multigaussFitDistribution,~,N] = multigauss(DipEvoFcn,K,r,2);

err(1) = any(abs(FitDistribution - multigaussFitDistribution)>1e-9);
err = any(err);

maxerr = max(abs(FitDistribution - multigaussFitDistribution));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(t,K*FitDistribution,'b')
   plot(t,K*multigaussFitDistribution,'r')
   subplot(122)
   hold on
   plot(r,FitDistribution,'b')
   plot(r,multigaussFitDistribution,'r')
   
end

end