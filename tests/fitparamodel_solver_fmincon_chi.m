function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [2 0.3];
P = rd_onegaussian(r,[2,0.3]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

[~,FitP] = fitparamodel(DipEvoFcn,@rd_onegaussian,r,K,0.6*InputParam,'Solver','fmincon','costmodel','chisquare');
err(1) = any(abs(FitP - P)>3e-7);
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