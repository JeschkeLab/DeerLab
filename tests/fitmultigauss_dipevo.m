function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [4 0.2 4 1 3 0.4 0.4 0.4];
P = rd_threegaussian(r,InputParam);


K = dipolarkernel(t,r);
DipEvoFcn = K*P;

[FitP,FitParam] = fitmultigauss(DipEvoFcn,K,r,5,'aicc');
err = any(abs(FitP - P)>7e-2);


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
end

end