function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [2.5 0.5 4 0.5 3 0.2 0.3 0.4];
P = rd_threegaussian(r,InputParam);

K = dipolarkernel(t,r);
B = exp(-0.15*t)';

lam = 0.5;
V = dipolarsignal(t,r,P,'moddepth',lam,'background',B);


KB = dipolarkernel(t,r,lam,B);

InitialGuess = [2 0.3 4 0.1 1 0.2 0.1 0.5];
[~,FitP] = fitparamodel(V,@rd_threegaussian,r,KB,InitialGuess);
err(1) = any(abs(FitP - P)>3e-4);
err = any(err);
maxerr = max(abs(FitP - P));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(t,V,'b')
   plot(t,KB*FitP,'r')
   subplot(122)
   hold on
   plot(r,P,'b')
   plot(r,FitP,'r')
   
end

end