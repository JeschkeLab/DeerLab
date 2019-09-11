function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.3 4 0.3 5 0.3 0.3 0.2];
P = rd_threegaussian(r,InputParam);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

[~,FitP] = fitparamodel(DipEvoFcn,@rd_threegaussian,r,K,'solver','lsqnonlin');
err(1) = any(abs(FitP - P)>1e-1);
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
   
end

end