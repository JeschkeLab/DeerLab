function [err,data,maxerr] = test(opt,oldata)

rng(2)
Dimension = 300;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.3]);

K = dipolarkernel(t,r);
DipEvoFcn = dipolarsignal(t,r,P,'noiselevel',0.01);

[~,FitP] = fitparamodel(DipEvoFcn,@rd_onegaussian,r,K,'costmodel','chisquare');
err(1) = any(abs(FitP - P)>8e-2);
err = any(err);

maxerr = max(abs(FitP - P));
data = [];

if opt.Display
     figure(1),clf
   subplot(121)
   hold on
   plot(t,DipEvoFcn,'b')
   plot(t,K*FitP,'r')
   subplot(122)
   hold on
   plot(r,P,'b')
   plot(r,FitP,'r')
   legend('truth','fit')
end

end