function [err,data,maxerr] = test(opt,oldata)


Dimension = 100;
t = linspace(0,5,Dimension);
r = time2dist(t);
InputParam = [4 0.2 4 1 3 0.4 0.4 0.4];
P = rd_threegaussian(r,InputParam);


B = td_exp(t,0.55);
V = dipolarsignal(t,r,P,'moddepth',0.5,'background',B);

[FitP,fitp] = fitmultigauss(V,t,r,5,'aicc','background',@td_exp);
err = any(abs(FitP - P)>8e-1);


maxerr = max(abs(FitP - P));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(t,V,'b')
   K = dipolarkernel(t,r,fitp(end-1),td_exp(t,fitp(end)));
   plot(t,K*FitP,'r')
   subplot(122)
   hold on
   plot(r,P,'b')
   plot(r,FitP,'r')
end

end